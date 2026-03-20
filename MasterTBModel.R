# TB Screening CEA — Markov cohort model (16 states, 43 strategies, lifetime horizon)
# MSc Modelling for Global Health, University of Oxford

library(ggplot2)    # plotting
library(patchwork)  # multi-panel figure assembly (must load before ggplot2 + overrides)
library(tidyr)      # data reshaping
library(dplyr)      # data manipulation
library(tibble)     # modern data frames
library(readxl)     # reading Excel initial conditions
library(ggrepel)    # non-overlapping text labels on plots
library(parallel)   # mclapply — PSA parallelisation

# Suppress scientific notation globally (never use e.g. 1e+05 in plots or output)
options(scipen = 999)

# -----------------------------------------------------------------------------
# Load all model parameters from config.csv
# Each row in config.csv specifies a parameter value, its probability
# distribution (for PSA), and its source reference.
# -----------------------------------------------------------------------------
config <- read.csv("config.csv", stringsAsFactors = FALSE)
config_list    <- setNames(config$value,                    config$parameter)
config_dist    <- setNames(config$distribution,             config$parameter)
config_dist_p1 <- setNames(as.numeric(config$dist_param1), config$parameter)
config_dist_p2 <- setNames(as.numeric(config$dist_param2), config$parameter)

# -----------------------------------------------------------------------------
# Model structural parameters
# Time horizon: 55 years (660 monthly cycles), entry age ~25 to age ~80 (lifetime; NICE reference case)
# Cohort size: 100,000 hypothetical migrants
# Discount rate: 3.5% per annum, applied to both costs and effects (NICE TA)
# -----------------------------------------------------------------------------
n_t      <- as.integer(config_list["n_t"])
n_c      <- as.integer(config_list["n_c"])
discount <- as.numeric(config_list["discount"])
discount_m   <- (1 + discount)^(1/12) - 1          # monthly compound discount rate
disc_factors <- 1 / (1 + discount_m)^seq_len(n_t)  # pre-computed 660-element discount vector

# ONS National Life Tables 2020-2022 (England & Wales, persons, both sexes).
# All-cause mortality rates per 100,000/year, by 5-year age band.
# Cohort entry age ~25; band 1 = 25-29, ..., band 11 = 75-79.
# 11 bands x 60 cycles (5yr) = 660 cycles = 55yr lifetime horizon.
# Converted to monthly probability: rate / 100,000 / 12.
# Source: ONS National Life Tables 2020-2022, Table 3 (England & Wales, persons).
# Values in config.csv (parameters ons_mort_25_29 through ons_mort_75_79).
ons_bg_mortality_monthly <- as.numeric(config_list[
  c("ons_mort_25_29","ons_mort_30_34","ons_mort_35_39","ons_mort_40_44",
    "ons_mort_45_49","ons_mort_50_54","ons_mort_55_59","ons_mort_60_64",
    "ons_mort_65_69","ons_mort_70_74","ons_mort_75_79")
]) / 1e5 / 12

# 16 health states representing the TB disease and care cascade.
# States are grouped into: uninfected, latent TB (LTBI), active TB, and dead.
# Within LTBI and active TB, states reflect diagnostic and treatment status.
v_state_names <- c(
  "Uninfected",
  "LatentUndiagnosed", "LatentDiagnosed",
  "ActiveUndiagnosed", "ActiveDiagnosed",
  "LatentTreated", "LatentNotreated",
  "ActiveTreated", "ActiveNotreated",
  "LatentCompleted", "LatentDiscontinued", "LatentLtfu",
  "ActiveCompleted", "ActiveDiscontinued", "ActiveLtfu",
  "Dead"
)

# Extract the subset of parameters used directly inside the model() function
# (state costs, QALY weights, transition probabilities). Other parameters
# such as initial conditions and effective detection rates are read separately.
param_names  <- names(config_list)
model_params <- param_names[grepl("^(Cost_|qaly_|p_)", param_names)]
paramsData   <- as.list(as.numeric(config_list[model_params]))
names(paramsData) <- model_params

# -----------------------------------------------------------------------------
# Core Markov model function
#
# Implements a monthly-cycle Markov cohort model with 16 health states.
# The model accepts an optional initial state distribution (init_dist) so
# that strategy-specific starting conditions (from the diagnostic decision
# tree) can be passed in. If no init_dist is provided, the model uses the
# background prevalence values from config.csv (no-screening baseline).
#
# Returns: state membership trace, discounted costs, discounted QALYs,
#          transition matrix, and per-state payoff vectors.
# -----------------------------------------------------------------------------
model <- function(.params, init_dist = NULL, n_t_override = NULL) {
  # Allow extended-horizon scenario analyses to shadow the global n_t without
  # modifying global state. If n_t_override is supplied, it is used for all
  # cycle counts in this call; otherwise the global n_t is used unchanged.
  if (!is.null(n_t_override)) n_t <- n_t_override
  with(.params, {
    n_s <- length(v_state_names)

    # Construct the 16x16 monthly transition probability matrix.
    # Rows = origin state, columns = destination state.
    # Off-diagonal entries are filled from config.csv parameters;
    # diagonal entries are set so each row sums to exactly 1.
    m_p <- matrix(0, nrow = n_s, ncol = n_s, dimnames = list(from = v_state_names, to = v_state_names))

    # -------------------------------------------------------------------------
    # Complete list of allowed state transitions (42 connections total).
    # All other cell entries remain zero (transitions not listed are forbidden).
    #
    # ENTRY (initial allocation at t = 0, handled outside model() via init_dist):
    #   U  -> U          uninfected remain uninfected
    #   L  -> Lu         latent TB missed by screening (proportion = 1 - eff_ltbi)
    #   L  -> Ld         latent TB detected by screening (proportion = eff_ltbi)
    #   I  -> Iu         active TB missed by screening (proportion = 1 - eff_active)
    #   I  -> Id         active TB detected by screening (proportion = eff_active)
    #
    # UNINFECTED (U):
    #   U  -> Lu         new LTBI acquisition (background transmission rate)
    #   U  -> Dead       background mortality (ONS)
    #
    # LATENT UNDIAGNOSED (Lu):
    #   Lu -> Ld         opportunistic/background screening detection
    #   Lu -> Iu         progression from latent to active TB (undiagnosed; observed incidence)
    #   Lu -> Dead       background mortality
    #
    # LATENT DIAGNOSED (Ld):
    #   Ld -> LTt        starts preventive treatment
    #   Ld -> LTn        declines / not offered treatment
    #   Ld -> Iu         reactivates despite diagnosis (before treatment starts)
    #   Ld -> Dead       background mortality
    #
    # LATENT ON TREATMENT (LTt):
    #   LTt -> LTOc      completes treatment course
    #   LTt -> LTOd      discontinues (adverse effects, ~9% UK)
    #   LTt -> LTOltfu   lost to follow-up during treatment
    #   LTt -> Dead      background + hepatotoxicity mortality
    #
    # LATENT NOT TREATED (LTn):
    #   LTn -> Iu        progression from latent to active TB (elevated risk, ~1%/year)
    #   LTn -> Dead      background mortality
    #
    # LATENT COMPLETED (LTOc):
    #   LTOc -> U        protected / effectively cured (treatment 60-90% effective)
    #   LTOc -> Dead     background mortality
    #
    # LATENT DISCONTINUED (LTOd):
    #   LTOd -> Iu       progression from latent to active TB (partial protection lost)
    #   LTOd -> Ld       re-engages with care
    #   LTOd -> Dead     background mortality
    #
    # LATENT LTFU (LTOltfu):
    #   LTOltfu -> Iu    progression from latent to active TB (no treatment protection)
    #   LTOltfu -> Ld    re-engages with care
    #   LTOltfu -> Dead  background mortality
    #
    # ACTIVE UNDIAGNOSED (Iu):
    #   Iu -> Id         diagnosed (spontaneous presentation / contact tracing)
    #   Iu -> Dead       elevated TB mortality (untreated ~50% 5-year mortality)
    #
    # ACTIVE DIAGNOSED (Id):
    #   Id -> ITt        starts treatment (2HRZE/4HR)
    #   Id -> ITn        declines / delays treatment
    #   Id -> Dead       pre-treatment mortality
    #
    # ACTIVE ON TREATMENT (ITt):
    #   ITt -> ITOc      completes treatment (84.4% UK 2023)
    #   ITt -> ITOd      discontinues treatment (~2% UK)
    #   ITt -> ITOltfu   lost to follow-up during treatment
    #   ITt -> Dead      treatment-phase mortality
    #
    # ACTIVE NOT TREATED (ITn):
    #   ITn -> Id        returns to care
    #   ITn -> Dead      high mortality without treatment
    #
    # ACTIVE COMPLETED (ITOc):
    #   ITOc -> U        cured (88% treatment success UK; post-TB sequelae captured via utility)
    #   ITOc -> Dead     post-treatment mortality (near background)
    #
    # ACTIVE DISCONTINUED (ITOd):
    #   ITOd -> Iu       returns to infectious / undiagnosed pool
    #   ITOd -> Dead     elevated mortality
    #
    # ACTIVE LTFU (ITOltfu):
    #   ITOltfu -> Iu    returns to infectious pool
    #   ITOltfu -> Dead  elevated mortality
    #
    # DEAD:
    #   Dead -> Dead     absorbing state
    # -------------------------------------------------------------------------

    # --- Disease progression and care cascade transitions ---
    m_p["Uninfected", "LatentUndiagnosed"] <- p_Uninfected_LatentUndiagnosed
    m_p["LatentUndiagnosed", "LatentDiagnosed"] <- p_LatentUndiagnosed_LatentDiagnosed
    m_p["LatentUndiagnosed", "ActiveUndiagnosed"] <- p_LatentUndiagnosed_ActiveUndiagnosed
    m_p["LatentDiagnosed", "ActiveUndiagnosed"] <- p_LatentDiagnosed_ActiveUndiagnosed
    m_p["ActiveUndiagnosed", "ActiveDiagnosed"] <- p_ActiveUndiagnosed_ActiveDiagnosed
    m_p["LatentDiagnosed", "LatentTreated"] <- p_LatentDiagnosed_LatentTreated
    m_p["LatentDiagnosed", "LatentNotreated"] <- p_LatentDiagnosed_LatentNotreated
    m_p["ActiveDiagnosed", "ActiveTreated"] <- p_ActiveDiagnosed_ActiveTreated
    m_p["ActiveDiagnosed", "ActiveNotreated"] <- p_ActiveDiagnosed_ActiveNotreated
    m_p["LatentTreated", "LatentCompleted"] <- p_LatentTreated_LatentCompleted
    m_p["LatentTreated", "LatentDiscontinued"] <- p_LatentTreated_LatentDiscontinued
    m_p["LatentTreated", "LatentLtfu"] <- p_LatentTreated_LatentLtfu
    m_p["ActiveTreated", "ActiveCompleted"] <- p_ActiveTreated_ActiveCompleted
    m_p["ActiveTreated", "ActiveDiscontinued"] <- p_ActiveTreated_ActiveDiscontinued
    m_p["ActiveTreated", "ActiveLtfu"] <- p_ActiveTreated_ActiveLtfu
    m_p["LatentNotreated", "ActiveUndiagnosed"] <- p_LatentNotreated_ActiveUndiagnosed
    m_p["ActiveDiscontinued", "ActiveUndiagnosed"] <- p_ActiveDiscontinued_ActiveUndiagnosed
    m_p["ActiveLtfu", "ActiveUndiagnosed"] <- p_ActiveLtfu_ActiveUndiagnosed
    m_p["LatentCompleted", "Uninfected"]        <- p_LatentCompleted_Uninfected
    # Residual reactivation after LTBI treatment: 14% of untreated rate (Berrocal-Almanza 2022: HR 0.14)
    m_p["LatentCompleted", "ActiveUndiagnosed"] <- p_LatentCompleted_ActiveUndiagnosed
    m_p["ActiveCompleted", "Uninfected"] <- p_ActiveCompleted_Uninfected

    # --- Progression from latent to active TB (observed incidence) and return-to-care transitions ---
    m_p["LatentDiscontinued", "ActiveUndiagnosed"] <- p_LatentDiscontinued_ActiveUndiagnosed
    m_p["LatentLtfu", "ActiveUndiagnosed"]         <- p_LatentLtfu_ActiveUndiagnosed
    m_p["ActiveNotreated", "ActiveDiagnosed"]      <- p_ActiveNotreated_ActiveDiagnosed

    # Re-engagement with care after LTBI treatment discontinuation or LTFU
    if ("p_LatentDiscontinued_LatentDiagnosed" %in% names(.params))
      m_p["LatentDiscontinued", "LatentDiagnosed"] <- .params[["p_LatentDiscontinued_LatentDiagnosed"]]
    if ("p_LatentLtfu_LatentDiagnosed" %in% names(.params))
      m_p["LatentLtfu", "LatentDiagnosed"] <- .params[["p_LatentLtfu_LatentDiagnosed"]]

    # --- Mortality transitions ---
    # State-specific monthly probabilities of death are read from config.csv.
    # Active TB states have elevated mortality; LTBI states use background rate.
    for(st in v_state_names[v_state_names != "Dead"]) {
      dead_par_name <- paste0("p_", st, "_Dead")
      if (dead_par_name %in% names(.params)) {
        m_p[st, "Dead"] <- .params[[dead_par_name]]
      }
    }
    # Dead is an absorbing state — once entered, no further transitions occur.
    m_p["Dead", "Dead"] <- 1

    # Diagonal = 1 − row sum (each row must sum to 1).
    off <- rowSums(m_p)
    for(i in seq_len(n_s)) {
      st <- v_state_names[i]
      if(st != "Dead") diag(m_p)[i] <- 1 - off[i]
    }

    # Numerical safeguard: during PSA, sampled parameter combinations can
    # occasionally cause row sums to exceed 1 (diagonal goes negative).
    # In this case, outgoing probabilities are rescaled proportionally so
    # the row still sums to 1. This is standard practice in PSA for Markov models.
    for(i in seq_len(n_s)) {
      st <- v_state_names[i]
      if(st == "Dead") next
      if(diag(m_p)[i] < 0) {
        off_diag_sum <- off[i]
        for(j in seq_len(n_s)) {
          if(i != j) m_p[i, j] <- m_p[i, j] / off_diag_sum
        }
        diag(m_p)[i] <- 0
      }
    }

    # --- Run the Markov cohort through n_t monthly cycles ---
    state_membership <- matrix(0, nrow = n_t, ncol = n_s, dimnames = list(1:n_t, v_state_names))

    # Populate cycle 1 from the strategy-specific initial distribution (if
    # provided) or from background prevalence values in config.csv (used for
    # the no-screening base case analysis).
    if (!is.null(init_dist)) {
      state_membership[1, ] <- init_dist
    } else {
      init_Uninfected        <- as.numeric(config_list["init_Uninfected"])
      init_LatentUndiagnosed <- as.numeric(config_list["init_LatentUndiagnosed"])
      init_ActiveUndiagnosed <- as.numeric(config_list["init_ActiveUndiagnosed"])
      state_membership[1, ]                        <- rep(0, n_s)
      state_membership[1, "Uninfected"]            <- n_c * init_Uninfected
      state_membership[1, "LatentUndiagnosed"]     <- n_c * init_LatentUndiagnosed
      state_membership[1, "ActiveUndiagnosed"]     <- n_c * init_ActiveUndiagnosed
    }

    # --- Age-varying background mortality (ONS 5-year age bands) ---
    # The cohort enters at age ~25. Background mortality increases with age;
    # the transition probability for background-mortality states is updated
    # every 60 cycles (5 years) to the relevant ONS age-band rate.
    # Active TB states keep their fixed TB-specific mortality (dominates background).
    # LatentTreated uses background + hepatotoxicity excess from NICE NG33.
    bg_mort_base     <- .params[["p_Uninfected_Dead"]]
    hepatotox_excess <- max(0, .params[["p_LatentTreated_Dead"]] - bg_mort_base)
    # ActiveCompleted is included so post-treatment mortality rises with age (ONS bands).
    # Omitting it caused a fixed 0.0002/mo rate that falls BELOW background at ages >65.
    bg_states <- c("Uninfected", "LatentUndiagnosed", "LatentDiagnosed",
                   "LatentNotreated", "LatentCompleted", "LatentDiscontinued", "LatentLtfu",
                   "ActiveCompleted")

    # Initialise matrix to age-band 1 (ages 25-29) before entering the loop.
    current_band <- 1L
    new_bg <- ons_bg_mortality_monthly[current_band]
    for (st in bg_states) {
      delta <- m_p[st, "Dead"] - new_bg
      m_p[st, "Dead"] <- new_bg
      m_p[st, st]     <- m_p[st, st] + delta
    }
    lt_delta <- m_p["LatentTreated", "Dead"] - (new_bg + hepatotox_excess)
    m_p["LatentTreated", "Dead"]          <- new_bg + hepatotox_excess
    m_p["LatentTreated", "LatentTreated"] <- m_p["LatentTreated", "LatentTreated"] + lt_delta

    # Two-phase reactivation: all latent→active transitions scale by
    # (p_react_phase2 / phase1_rate) at cycle 61 (year 5).
    # To disable: set p_react_phase2 = p_LatentUndiagnosed_ActiveUndiagnosed in config.csv.
    react_phase1         <- .params[["p_LatentUndiagnosed_ActiveUndiagnosed"]]
    p_react_phase2       <- .params[["p_react_phase2"]]
    react_phase_switched <- FALSE
    react_latent_states  <- c("LatentUndiagnosed", "LatentDiagnosed", "LatentNotreated",
                               "LatentCompleted",   "LatentDiscontinued", "LatentLtfu")

    # Matrix multiplication advances the cohort one cycle at a time.
    # At each 5-year band boundary (cycles 61, 121, ..., 601) the background
    # mortality rows are updated to the next ONS age-band rate (11 bands total,
    # covering ages 25-79; band 11 applies for ages 75-80 at end of horizon).
    for(t in 2:n_t) {
      new_band <- min(11L, as.integer(floor((t - 1) / 60)) + 1L)
      if (new_band != current_band) {
        current_band <- new_band
        new_bg <- ons_bg_mortality_monthly[current_band]
        for (st in bg_states) {
          delta <- m_p[st, "Dead"] - new_bg
          m_p[st, "Dead"] <- new_bg
          m_p[st, st]     <- m_p[st, st] + delta
        }
        lt_delta <- m_p["LatentTreated", "Dead"] - (new_bg + hepatotox_excess)
        m_p["LatentTreated", "Dead"]          <- new_bg + hepatotox_excess
        m_p["LatentTreated", "LatentTreated"] <- m_p["LatentTreated", "LatentTreated"] + lt_delta
      }
      # Switch reactivation to phase 2 rate at year 5 (cycle 61)
      if (!react_phase_switched && t == 61L) {
        react_phase_switched <- TRUE
        react_scale <- p_react_phase2 / react_phase1
        for (st in react_latent_states) {
          old_r            <- m_p[st, "ActiveUndiagnosed"]
          new_r            <- old_r * react_scale
          m_p[st, "ActiveUndiagnosed"] <- new_r
          m_p[st, st]      <- m_p[st, st] + (old_r - new_r)
        }
      }
      state_membership[t, ] <- as.numeric(state_membership[t-1, ] %*% m_p)
      if (abs(sum(state_membership[t, ]) - n_c) > 1.0)
        warning(sprintf("Cycle %d: cohort sum = %.1f (expected %.0f; check transition matrix row sums)", t, sum(state_membership[t, ]), n_c))
    }

    # --- Per-state monthly cost and QALY payoff vectors ---
    # Costs are in 2024 GBP. QALYs are divided by 12 to
    # convert annual utility weights to monthly values.
    m_payoffsCost <- matrix(0, nrow = n_s, ncol = 1, dimnames = list(v_state_names, "Costs"))
    m_payoffsQALY <- matrix(0, nrow = n_s, ncol = 1, dimnames = list(v_state_names, "QALYs"))
    m_payoffsCost[,"Costs"] <- c(
      Cost_Uninfected,
      Cost_LatentUndiagnosed,Cost_LatentDiagnosed,
      Cost_ActiveUndiagnosed,Cost_ActiveDiagnosed,
      Cost_LatentTreated,Cost_LatentNotreated,
      Cost_ActiveTreated,Cost_ActiveNotreated,
      Cost_LatentCompleted,Cost_LatentDiscontinued,Cost_LatentLtfu,
      Cost_ActiveCompleted,Cost_ActiveDiscontinued,Cost_ActiveLtfu,
      Cost_Dead
    )
    m_payoffsQALY[,"QALYs"] <- c(
      qaly_Uninfected,
      qaly_LatentUndiagnosed,qaly_LatentDiagnosed,
      qaly_ActiveUndiagnosed,qaly_ActiveDiagnosed,
      qaly_LatentTreated,qaly_LatentNotreated,
      qaly_ActiveTreated,qaly_ActiveNotreated,
      qaly_LatentCompleted,qaly_LatentDiscontinued,qaly_LatentLtfu,
      qaly_ActiveCompleted,qaly_ActiveDiscontinued,qaly_ActiveLtfu,
      qaly_Dead
    ) / 12

    # Apply discounting at 3.5% per annum (converted to monthly compound rate)
    # per NICE reference case. Future costs and QALYs are worth less than
    # immediate ones; dividing by (1 + r)^t adjusts for this.
    payoff_traceCost <- state_membership %*% m_payoffsCost
    payoff_traceQaly <- state_membership %*% m_payoffsQALY

    # Vectorised discounting: multiply by pre-computed discount vector.
    # If n_t was overridden (scenario analysis), compute locally; otherwise use
    # the global disc_factors (660-element vector, computed once at startup).
    dv <- if (n_t == length(disc_factors)) disc_factors else 1/(1+discount_m)^seq_len(n_t)
    payoff_traceCost_d <- payoff_traceCost * dv
    payoff_traceQaly_d <- payoff_traceQaly * dv

    # Half-cycle correction: assumes transitions occur on average mid-cycle
    # rather than at the start. Subtracts half the first and last cycle values
    # to correct for the overcount introduced by the standard Markov assumption.
    total_cost_hcc <- sum(payoff_traceCost_d) -
      0.5 * payoff_traceCost_d[1,1] - 0.5 * payoff_traceCost_d[n_t,1]
    total_qaly_hcc <- sum(payoff_traceQaly_d) -
      0.5 * payoff_traceQaly_d[1,1] - 0.5 * payoff_traceQaly_d[n_t,1]

    list(
      state_membership = state_membership,
      payoff_trace_perCycle_cost = payoff_traceCost_d,
      payoff_trace_perCycle_qaly = payoff_traceQaly_d,
      total_cost = total_cost_hcc,
      total_qaly = total_qaly_hcc,
      m_p = m_p,
      m_payoffsCost = m_payoffsCost,
      m_payoffsQALY = m_payoffsQALY
    )
  })
}

# -----------------------------------------------------------------------------
# PSA parameter sampling
#
# Draws one random sample of all model parameters from the probability
# distributions specified in config.csv. Beta distributions are used for
# probabilities and utilities (bounded 0-1); gamma distributions are used
# for costs (bounded >0, right-skewed). Fixed parameters are not varied.
# This function is called once per PSA simulation inside run_psa().
# -----------------------------------------------------------------------------
sample_params <- function() {
  sampled <- paramsData  # start from base case values

  # Vectorised sampling: draw all beta and gamma parameters in two batch calls.
  # R's rbeta/rgamma recycle shape vectors, sampling one value per parameter set.
  beta_names  <- names(sampled)[config_dist[names(sampled)] == "beta"  & !is.na(config_dist[names(sampled)])]
  gamma_names <- names(sampled)[config_dist[names(sampled)] == "gamma" & !is.na(config_dist[names(sampled)])]

  if (length(beta_names) > 0) {
    sampled[beta_names] <- as.list(rbeta(length(beta_names),
                                         shape1 = config_dist_p1[beta_names],
                                         shape2 = config_dist_p2[beta_names]))
  }
  if (length(gamma_names) > 0) {
    sampled[gamma_names] <- as.list(rgamma(length(gamma_names),
                                            shape = config_dist_p1[gamma_names],
                                            scale = config_dist_p2[gamma_names]))
  }

  return(sampled)
}

#------------------------------------------------------------------------------#
# Diagnostic strategy initial conditions
#
# Active TB allocation: TP and FN taken directly from the Zenner et al. 2025 ERJ
#   decision tree (Initial_conditions.xlsx, 100k cohort). These account for the
#   full pathway (symptom screen → test → culture confirmation).
#   TP → ActiveDiagnosed | FN → ActiveUndiagnosed
#   FP + TN → Uninfected (UK assumption: culture confirmation before treatment)
#
# LTBI allocation: taken directly from Zenner et al. 2025 ERJ decision tree (Excel col 5).
#   LatentUndiagnosed = 17,800 (Berrocal-Almanza 2022: 17.8% per 100k) − Excel detected.
#   Excel LTBI column used directly; consistent with decision tree design.
#   The decision tree embeds its own lower LTBI prevalence assumption (~0.1–0.15%),
#   which differs from the Berrocal-Almanza 2022 17.8% background prevalence used elsewhere.
#
# Upfront diagnostic costs: Tot_costs from Excel, per-person = Tot_costs / 100,000.
#------------------------------------------------------------------------------#

# Background prevalence at entry — read from config.csv (UKHSA 2024)
prev_uninfected <- as.numeric(config_list["init_Uninfected"])
prev_ltbi       <- as.numeric(config_list["init_LatentUndiagnosed"])
prev_active     <- as.numeric(config_list["init_ActiveUndiagnosed"])

# IGRA LTBI sensitivities — read from config.csv (Zenner et al. 2025 Eur Respir J Table 1)
igra_sens_qft  <- as.numeric(config_list["igra_ltbi_sens_qft"])   # QFT-TB Gold Plus: 0.83
igra_sens_tspt <- as.numeric(config_list["igra_ltbi_sens_tspt"])  # T-SPOT.TB: 0.88

# IGRA LTBI specificities — read from config.csv (Pai et al. 2008, BCG-vaccinated populations)
# QFT: 0.96 (95% CI 0.94-0.98); T-SPOT: 0.93 (95% CI 0.86-1.00)
igra_spec_qft  <- as.numeric(config_list["igra_spec_qft"])        # QFT specificity:   0.96
igra_spec_tspt <- as.numeric(config_list["igra_spec_tspt"])       # T-SPOT specificity: 0.93

# eff_active parameter names (for PSA sampling — all have beta distributions in config.csv)
eff_active_param_names <- c("eff_active_CXR_only", "eff_active_QFT",
                              "eff_active_TSPOT",   "eff_active_TST",
                              "eff_active_CXR_QFT", "eff_active_TST_QFT",
                              "eff_active_CXR_Xpert")
eff_active_base <- setNames(
  as.numeric(config_list[eff_active_param_names]),
  eff_active_param_names
)
eff_active_dist_p1 <- setNames(
  as.numeric(config_dist_p1[eff_active_param_names]),
  eff_active_param_names
)
eff_active_dist_p2 <- setNames(
  as.numeric(config_dist_p2[eff_active_param_names]),
  eff_active_param_names
)

# Read decision tree outputs from Excel (active TB TP/FN and diagnostic costs)
# Columns (after skipping 2 header rows): 1=strategy, 4=TP, 5=LTBI, 6=Tot_costs, 9=FN
excel_ic <- read_excel("Supplementary files/S2. Screening decision tree - Initial conditions for Markov model.xlsx", col_names = FALSE, skip = 2)
excel_denom <- 100000  # Excel cohort denominator (TP+FN+TN+FP = 100,000)

# Column index validation — guards against silent errors if Excel columns shift
# Expected layout (after skip=2): col1=strategy, col4=TP, col5=LTBI, col6=Tot_costs, col9=FN
# Col 4 (TP): active TB true positives per 100k — expect 100–1,000
# Col 6 (Tot_costs): programme costs per 100k — expect >100,000 (typically 400k–19M)
# Col 9 (FN): active TB false negatives per 100k — expect 0–900
stopifnot(
  "Excel col 1 (strategy names): first entry should be a string" =
    is.character(excel_ic[[1]][1]),
  "Excel col 4 (TP): values should be 100-1000 (active TB true positives per 100k)" =
    all(na.omit(as.numeric(excel_ic[[4]])) > 50 &
        na.omit(as.numeric(excel_ic[[4]])) < 1500),
  "Excel col 6 (Tot_costs): values should be >100,000 (programme costs per 100k)" =
    median(na.omit(as.numeric(excel_ic[[6]]))) > 100000,
  "Excel col 9 (FN): values should be 0-900 (active TB false negatives per 100k)" =
    all(na.omit(as.numeric(excel_ic[[9]])) >= 0 &
        na.omit(as.numeric(excel_ic[[9]])) < 1000)
)
cat("Excel column validation passed.\n")

get_ic <- function(row_name) {
  r <- excel_ic[excel_ic[[1]] == row_name, ]
  list(
    TP         = as.numeric(r[[4]]),
    LTBI       = as.numeric(r[[5]]),
    Tot_costs  = as.numeric(r[[6]]),
    FN         = as.numeric(r[[9]])
  )
}


# Project-wide colour palette — pink-purple scheme throughout
project_pal <- colorRampPalette(c(
  "#2d0a2e", "#5a1a7a", "#7a1a4a", "#a03060",
  "#c84b6a", "#e0788a", "#e8a0b0", "#f5d0dc"
))

# Build named colour vector for any set of strategy names
build_strategy_colours <- function(nms) {
  n <- length(nms)
  cols <- project_pal(n)
  setNames(cols, sort(nms))
}

# Clean Excel row names -> human-readable labels for plots
clean_excel_name <- function(nm) {
  known <- c(
    "anycough_cxr(TB)"           = "Cough+CXR (TB sx)",
    "anysx_cxr(any)"             = "Symptom screen+CXR",
    "cough/CXR (parallel)_xpert" = "CXR+Xpert",
    "cough/CXR (parallel)_ultra" = "CXR+Ultra",
    "qft_xpert"                  = "QFT-GIT+Xpert",
    "tspt_xpert"                 = "T-SPOT.TB+Xpert",
    "tst_cxr(TB)"                = "CXR+TST",
    "qft_cxr(TB)"                = "CXR+QFT-GIT",
    "parallel allsx qft_xpert"   = "Parallel Sx+QFT (Xpert)",
    "parallel allsx qft_ultra"   = "Parallel Sx+QFT (Ultra)",
    "parallel cough qft_xpert"   = "Parallel Cough+QFT (Xpert)",
    "parallel cough qft_ultra"   = "Parallel Cough+QFT (Ultra)",
    "parallel allsx tspt_xpert"  = "Parallel Sx+T-SPOT (Xpert)",
    "parallel allsx tspt_ultra"  = "Parallel Sx+T-SPOT (Ultra)",
    "parallel cough tspt_xpert"  = "Parallel Cough+T-SPOT (Xpert)",
    "parallel cough tspt_ultra"  = "Parallel Cough+T-SPOT (Ultra)"
  )
  if (nm %in% names(known)) return(known[[nm]])
  out <- nm
  out <- gsub("cough/CXR \\(parallel\\)", "CXR", out)
  out <- gsub("^anycough", "Cough", out)
  out <- gsub("^anysx",    "Sx screen", out)
  out <- gsub("_", " + ", out)
  out <- gsub("\\btspt\\b|\\btspot\\b", "T-SPOT.TB", out, perl = TRUE)
  out <- gsub("\\bqft\\b",   "QFT-GIT",  out, perl = TRUE, ignore.case = TRUE)
  out <- gsub("\\btst\\b",   "TST",      out, perl = TRUE, ignore.case = TRUE)
  out <- gsub("\\bcxr\\b",   "CXR",      out, perl = TRUE, ignore.case = TRUE)
  out <- gsub("\\bxpert\\b", "Xpert",    out, perl = TRUE, ignore.case = TRUE)
  out <- gsub("\\bultra\\b", "Ultra",    out, perl = TRUE, ignore.case = TRUE)
  out <- gsub("\\(TB\\)",    "(TB sx)",  out)
  out <- gsub("([[:alpha:]])\\(", "\\1 (", out)  # ensure space before (
  out <- gsub("\\s*\\+\\s*", "+", out)
  trimws(out)
}

# Assign test family by primary (first-line) test, not confirmatory test.
# Order matters: check IGRA/TST families before Xpert/Ultra so that
# e.g. "QFT-GIT+Xpert" is classified as QFT-GIT, not Molecular.
assign_family <- function(nm) {
  nm_lc <- tolower(nm)
  if (grepl("qft",         nm_lc))  return("QFT-GIT")
  if (grepl("tspot|tspt",  nm_lc))  return("T-SPOT.TB")
  if (grepl("^tst|_tst",   nm_lc))  return("TST Mantoux")
  if (grepl("xpert|ultra", nm_lc))  return("Molecular (Xpert/Ultra)")
  return("CXR / Symptom screen")
}

# eff_active key map: maps Excel row names to the eff_active_* config parameter
# that governs active TB detection uncertainty for that test pathway.
eff_active_key_map <- c(
  "anycough_cxr(any)"                       = "eff_active_CXR_only",
  "2 wk cough_cxr (any)"                    = "eff_active_CXR_only",
  "2 wk cough_cxr (TB)"                     = "eff_active_CXR_only",
  "anycough_cxr(TB)"                        = "eff_active_CXR_only",
  "anysx_cxr(any)"                          = "eff_active_CXR_only",
  "anysx_cxr(TB)"                           = "eff_active_CXR_only",
  "cough/CXR (parallel)_xpert"              = "eff_active_CXR_Xpert",
  "cough/CXR (parallel)_ultra"              = "eff_active_CXR_Xpert",
  "anycough_cxr(any)_xpert"                 = "eff_active_CXR_Xpert",
  "anycough_cxr(any)_ultra"                 = "eff_active_CXR_Xpert",
  "qft_cxr(any)"                            = "eff_active_QFT",
  "qft_cxr(TB)"                             = "eff_active_QFT",
  "qft_xpert"                               = "eff_active_QFT",
  "qft_ultra"                               = "eff_active_QFT",
  "qft_cxr(TB)_xpert"                       = "eff_active_CXR_QFT",
  "qft_cxr(TB)_ultra"                       = "eff_active_CXR_QFT",
  "qft_cough/CXR (parallel)_xpert"          = "eff_active_CXR_QFT",
  "qft_cough/CXR (parallel)_ultra"          = "eff_active_CXR_QFT",
  "tspt_cxr(any)"                           = "eff_active_TSPOT",
  "tspt_cxr(TB)"                            = "eff_active_TSPOT",
  "tspt_xpert"                              = "eff_active_TSPOT",
  "tspt_ultra"                              = "eff_active_TSPOT",
  "tspot_cxr(TB)_xpert"                     = "eff_active_TSPOT",
  "tspot_cxr(TB)_ultra"                     = "eff_active_TSPOT",
  "tspot_cough/CXR (parallel)_xpert"        = "eff_active_TSPOT",
  "tspot_cough/CXR (parallel)_ultra"        = "eff_active_TSPOT",
  "tst_cxr(any)"                            = "eff_active_TST",
  "tst_cxr(TB)"                             = "eff_active_TST",
  "tst_xpert"                               = "eff_active_TST",
  "tst_ultra"                               = "eff_active_TST",
  "tst_cxr(TB)_xpert"                       = "eff_active_TST_QFT",
  "tst_cxr(TB)_ultra"                       = "eff_active_TST_QFT",
  "tst_cough/CXR (parallel)_xpert"          = "eff_active_TST_QFT",
  "tst_cough/CXR (parallel)_ultra"          = "eff_active_TST_QFT",
  "parallel allsx qft_xpert"                = "eff_active_QFT",
  "parallel allsx qft_ultra"                = "eff_active_QFT",
  "parallel cough qft_xpert"                = "eff_active_QFT",
  "parallel cough qft_ultra"                = "eff_active_QFT",
  "parallel allsx tspt_xpert"               = "eff_active_TSPOT",
  "parallel allsx tspt_ultra"               = "eff_active_TSPOT",
  "parallel cough tspt_xpert"               = "eff_active_TSPOT",
  "parallel cough tspt_ultra"               = "eff_active_TSPOT"
)

create_diagnostic_strategies <- function() {
  make_init <- function() setNames(rep(0, length(v_state_names)), v_state_names)

  # Total LTBI and active TB in model cohort (from prevalence assumptions)
  n_ltbi   <- n_c * prev_ltbi    # Berrocal-Almanza 2022: 17.8% → 17,800 per 100k
  n_active <- n_c * prev_active  # 1% → 1,000 per 100k

  # Helper: build init vector from Excel TP/FN/LTBI columns directly
  make_from_excel <- function(ic) {
    init <- make_init()
    ltbi_detected               <- ic$LTBI * (n_c / excel_denom)
    init["ActiveDiagnosed"]     <- ic$TP * (n_c / excel_denom)
    init["ActiveUndiagnosed"]   <- ic$FN * (n_c / excel_denom)
    init["LatentDiagnosed"]     <- ltbi_detected
    init["LatentUndiagnosed"]   <- n_ltbi - ltbi_detected
    init["Uninfected"]          <- n_c * prev_uninfected
    init
  }

  # ---------------------------------------------------------------------------
  # All 34 feasible pathways from Zenner et al. 2025 ERJ decision tree
  # (rows with non-NA Tot_costs) plus No_screening reference = 35 total.
  # 8 parallel IGRA strategies with Option B LTBI override added below = 43 total.
  # ---------------------------------------------------------------------------

  # ---- Reference: No screening ----
  s <- list(name = "No screening", excel_name = NA_character_,
            init = make_init(), test_cost = 0)
  s$init["Uninfected"]        <- n_c * prev_uninfected
  s$init["LatentUndiagnosed"] <- n_c * prev_ltbi
  s$init["ActiveUndiagnosed"] <- n_c * prev_active
  strategies <- list(No_screening = s)
  strategies[["No_screening"]]$eff_active_key <- NA_character_
  strategies[["No_screening"]]$base_tp        <- 0
  strategies[["No_screening"]]$base_fn        <- n_active  # all active TB undetected

  # ---- All Excel rows with non-NA Tot_costs ----
  for (i in seq_len(nrow(excel_ic))) {
    row      <- excel_ic[i, ]
    excel_nm <- as.character(row[[1]])
    if (is.na(excel_nm) || excel_nm == "") next
    tc <- suppressWarnings(as.numeric(row[[6]]))
    if (is.na(tc)) next
    ic <- list(
      TP        = as.numeric(row[[4]]),
      LTBI      = as.numeric(row[[5]]),
      Tot_costs = tc,
      FN        = as.numeric(row[[9]])
    )
    key <- gsub("[^A-Za-z0-9]", "_", excel_nm)
    strategies[[key]] <- list(
      name       = clean_excel_name(excel_nm),
      excel_name = excel_nm,
      init       = make_from_excel(ic),
      test_cost  = ic$Tot_costs / excel_denom
    )
    strategies[[key]]$eff_active_key <- unname(eff_active_key_map[excel_nm])
    strategies[[key]]$base_tp        <- ic$TP * (n_c / excel_denom)
    strategies[[key]]$base_fn        <- ic$FN * (n_c / excel_denom)
  }
  # ---------------------------------------------------------------------------
  # OPTION B: 8 parallel IGRA strategies — LTBI recalculated from
  # population prevalence × pooled IGRA sensitivity (Zenner 2025 Table 1).
  # Excel LTBI values for these rows reflect sequential pathway yield and are
  # overridden here with universal-IGRA estimates.
  #
  # LTBI override: n_ltbi × igra_sensitivity
  #   QFT pooled sensitivity:     0.83 (Zenner 2025 Table 1)
  #   T-SPOT.TB sensitivity:      0.88 (Zenner 2025 Table 1)
  #
  # Cost estimation: existing IGRA+confirmatory base cost + IGRA testing programme cost
  #   allsx (any TB symptom) screen: cost_igra_programme_allsx per 100k (base £500,000; £5/pp)
  #   cough screen only:             cost_igra_programme_cough per 100k (base £300,000; £3/pp)
  #
  # Programme delivery cost values are read from paramsData (config.csv rows
  # cost_igra_programme_cough / cost_igra_programme_allsx) so they can be varied in PSA and DSA.
  # ---------------------------------------------------------------------------

  prog_cost_cough <- as.numeric(config_list["cost_igra_programme_cough"])  # £ per cohort of 100k
  prog_cost_allsx <- as.numeric(config_list["cost_igra_programme_allsx"])

  # Net LTBI detected = TP − FP
  # FP = (1 - specificity) × n_uninfected; source: Pai et al. 2008 BCG-vaccinated meta-analysis
  # QFT: 17,800 × 0.83 − 0.04 × 81,200 = 14,774 − 3,248 = 11,526/100k
  # T-SPOT: 17,800 × 0.88 − 0.07 × 81,200 = 15,664 − 5,684 = 9,980/100k
  igra_ltbi_qft  <- max(0, n_ltbi * igra_sens_qft  - (1 - igra_spec_qft)  * n_c * prev_uninfected)
  igra_ltbi_tspt <- max(0, n_ltbi * igra_sens_tspt - (1 - igra_spec_tspt) * n_c * prev_uninfected)

  # Read TP and FN for parallel IGRA strategies from Excel rows.
  # These rows (parallel allsx/cough qft/tspt xpert/ultra) have TP/FN populated
  # in the decision tree but Tot_costs = NA (hence excluded from the main loop).
  # If a row is missing or has NA values, falls back to Zenner 2025 Table 1 values.
  lookup_parallel_tp_fn <- function(nm, tp_fallback, fn_fallback) {
    r <- excel_ic[!is.na(excel_ic[[1]]) & excel_ic[[1]] == nm, ]
    if (nrow(r) == 0 || is.na(as.numeric(r[[4]])))  {
      warning(paste0("Excel row '", nm, "' not found or TP is NA — using Zenner 2025 Table 1 fallback: TP=", tp_fallback))
      return(list(tp = tp_fallback, fn = fn_fallback))
    }
    list(tp = as.numeric(r[[4]]), fn = as.numeric(r[[9]]))
  }

  make_parallel_igra_init <- function(tp, fn, ltbi_detected) {
    init <- make_init()
    init["ActiveDiagnosed"]   <- tp * (n_c / excel_denom)
    init["ActiveUndiagnosed"] <- fn * (n_c / excel_denom)
    init["LatentDiagnosed"]   <- ltbi_detected
    init["LatentUndiagnosed"] <- n_ltbi - ltbi_detected
    init["Uninfected"]        <- n_c * prev_uninfected
    init
  }

  qft_xpert_base  <- as.numeric(excel_ic[excel_ic[[1]] == "qft_xpert",  ][[6]])
  qft_ultra_base  <- as.numeric(excel_ic[excel_ic[[1]] == "qft_ultra",  ][[6]])
  tspt_xpert_base <- as.numeric(excel_ic[excel_ic[[1]] == "tspt_xpert", ][[6]])
  tspt_ultra_base <- as.numeric(excel_ic[excel_ic[[1]] == "tspt_ultra", ][[6]])

  parallel_igra_specs <- list(
    list(nm = "parallel allsx qft_xpert",  pf = lookup_parallel_tp_fn("parallel allsx qft_xpert",  836.6, 163.4), ltbi = igra_ltbi_qft,  cost = qft_xpert_base  + prog_cost_allsx, prog_cost = prog_cost_allsx, prog_cost_type = "allsx"),
    list(nm = "parallel allsx qft_ultra",  pf = lookup_parallel_tp_fn("parallel allsx qft_ultra",  874.6, 125.4), ltbi = igra_ltbi_qft,  cost = qft_ultra_base  + prog_cost_allsx, prog_cost = prog_cost_allsx, prog_cost_type = "allsx"),
    list(nm = "parallel cough qft_xpert",  pf = lookup_parallel_tp_fn("parallel cough qft_xpert",  791.7, 208.3), ltbi = igra_ltbi_qft,  cost = qft_xpert_base  + prog_cost_cough, prog_cost = prog_cost_cough, prog_cost_type = "cough"),
    list(nm = "parallel cough qft_ultra",  pf = lookup_parallel_tp_fn("parallel cough qft_ultra",  827.7, 172.3), ltbi = igra_ltbi_qft,  cost = qft_ultra_base  + prog_cost_cough, prog_cost = prog_cost_cough, prog_cost_type = "cough"),
    list(nm = "parallel allsx tspt_xpert", pf = lookup_parallel_tp_fn("parallel allsx tspt_xpert", 849.4, 150.6), ltbi = igra_ltbi_tspt, cost = tspt_xpert_base + prog_cost_allsx, prog_cost = prog_cost_allsx, prog_cost_type = "allsx"),
    list(nm = "parallel allsx tspt_ultra", pf = lookup_parallel_tp_fn("parallel allsx tspt_ultra", 888.0, 112.0), ltbi = igra_ltbi_tspt, cost = tspt_ultra_base + prog_cost_allsx, prog_cost = prog_cost_allsx, prog_cost_type = "allsx"),
    list(nm = "parallel cough tspt_xpert", pf = lookup_parallel_tp_fn("parallel cough tspt_xpert", 817.7, 182.3), ltbi = igra_ltbi_tspt, cost = tspt_xpert_base + prog_cost_cough, prog_cost = prog_cost_cough, prog_cost_type = "cough"),
    list(nm = "parallel cough tspt_ultra", pf = lookup_parallel_tp_fn("parallel cough tspt_ultra", 854.9, 145.1), ltbi = igra_ltbi_tspt, cost = tspt_ultra_base + prog_cost_cough, prog_cost = prog_cost_cough, prog_cost_type = "cough")
  )

  for (ps in parallel_igra_specs) {
    key <- gsub("[^A-Za-z0-9]", "_", ps$nm)
    # Determine eff_active_key from the strategy name:
    # QFT-based parallel strategies use eff_active_QFT; T-SPOT-based use eff_active_TSPOT.
    ps_eff_key <- if (grepl("tspt|tspot", ps$nm, ignore.case = TRUE)) "eff_active_TSPOT" else "eff_active_QFT"
    strategies[[key]] <- list(
      name            = clean_excel_name(ps$nm),
      excel_name      = ps$nm,
      init            = make_parallel_igra_init(ps$pf$tp, ps$pf$fn, ps$ltbi),
      test_cost       = ps$cost / excel_denom,
      prog_cost_total = ps$prog_cost,      # total programme delivery cost £ for cohort (used in DSA)
      prog_cost_type  = ps$prog_cost_type, # "cough" or "allsx"
      eff_active_key  = ps_eff_key,
      base_tp         = ps$pf$tp * (n_c / excel_denom),
      base_fn         = ps$pf$fn * (n_c / excel_denom)
    )
  }

  cat(sprintf("Parallel IGRA strategies added (Option B): LTBI detection = %.0f (QFT) / %.0f (T-SPOT) per 100k\n",
              igra_ltbi_qft, igra_ltbi_tspt))

  return(strategies)
}

# -----------------------------------------------------------------------------
# Run the Markov model for a single diagnostic strategy
#
# Passes the strategy-specific initial state distribution to model() and
# adds the one-time screening test cost to the total (test costs are incurred
# at entry and are not part of the recurring state-based cost structure).
# -----------------------------------------------------------------------------
run_strategy <- function(strategy, params = paramsData, n_t_override = NULL) {
  result <- model(params, init_dist = strategy$init, n_t_override = n_t_override)
  result$total_cost    <- result$total_cost + (strategy$test_cost * n_c)
  result$test_cost     <- strategy$test_cost * n_c
  result$strategy_name <- strategy$name
  return(result)
}

# -----------------------------------------------------------------------------
# Incremental cost-effectiveness ratio (ICER) calculation
#
# Performs a full sequential incremental analysis following NICE methods:
#   1. Pairwise ICERs vs No Screening (retained for reference).
#   2. Simple dominance: any strategy with higher cost AND lower QALYs than a
#      cheaper alternative is flagged "simply dominated".
#   3. Extended dominance: strategies whose sequential ICER exceeds that of the
#      next comparator on the efficient frontier are iteratively removed.
#   4. Sequential ICERs on the efficient frontier are reported.
#
# The NICE willingness-to-pay threshold is £25,000–£35,000 per QALY (updated Dec 2025, effective Apr 2026).
# -----------------------------------------------------------------------------
calculate_icer <- function(strategy_results) {
  df <- tibble(
    strategy = sapply(strategy_results, function(x) x$strategy_name),
    cost     = sapply(strategy_results, function(x) x$total_cost),
    qaly     = sapply(strategy_results, function(x) x$total_qaly),
    deaths   = sapply(strategy_results, function(x) x$state_membership[nrow(x$state_membership), "Dead"])
  )

  df <- df %>% arrange(cost)

  ref_idx <- which(df$strategy == "No screening")
  if (length(ref_idx) == 0) ref_idx <- 1

  # Pairwise ICERs vs No Screening (retained for backward compatibility)
  df <- df %>%
    mutate(
      inc_cost        = cost - cost[ref_idx],
      inc_qaly        = qaly - qaly[ref_idx],
      icer            = ifelse(inc_qaly != 0, inc_cost / inc_qaly, NA),
      cost_per_person = cost / n_c,
      qaly_per_person = qaly / n_c
    )

  # -------------------------------------------------------------------
  # Sequential incremental analysis with extended dominance detection
  # -------------------------------------------------------------------
  df$dominance      <- "non-dominated"
  df$sequential_icer <- NA_real_
  df$dominance[ref_idx] <- "ref"

  # Step 1: simple dominance — flag any strategy where a cheaper alternative
  #         achieves at least as many QALYs.
  for (i in seq_len(nrow(df))) {
    if (df$dominance[i] == "ref") next
    cheaper <- df[df$cost < df$cost[i], ]
    if (nrow(cheaper) > 0 && any(cheaper$qaly >= df$qaly[i])) {
      df$dominance[i] <- "simply dominated"
    }
  }

  # Step 2: extended dominance — iteratively remove strategies whose sequential
  #         ICER exceeds that of the next comparator (frontier not convex).
  changed <- TRUE
  while (changed) {
    changed  <- FALSE
    frontier <- df[df$dominance %in% c("ref", "non-dominated"), ] %>% arrange(cost)
    n_f      <- nrow(frontier)
    if (n_f < 3) break  # need at least 3 points to identify extended dominance
    for (i in 2:(n_f - 1)) {
      d_cost_i    <- frontier$cost[i]   - frontier$cost[i - 1]
      d_qaly_i    <- frontier$qaly[i]   - frontier$qaly[i - 1]
      d_cost_next <- frontier$cost[i + 1] - frontier$cost[i]
      d_qaly_next <- frontier$qaly[i + 1] - frontier$qaly[i]
      if (d_qaly_i <= 0 || d_qaly_next <= 0) next
      seq_icer_i    <- d_cost_i    / d_qaly_i
      seq_icer_next <- d_cost_next / d_qaly_next
      if (seq_icer_i > seq_icer_next) {
        idx <- which(df$strategy == frontier$strategy[i])
        df$dominance[idx] <- "extendedly dominated"
        changed <- TRUE
        break  # restart after each removal
      }
    }
  }

  # Step 3: assign sequential ICERs on the final efficient frontier.
  frontier_final <- df[df$dominance %in% c("ref", "non-dominated"), ] %>% arrange(cost)
  for (i in seq_len(nrow(frontier_final))) {
    idx <- which(df$strategy == frontier_final$strategy[i])
    if (i == 1) {
      df$sequential_icer[idx] <- NA_real_
    } else {
      d_cost <- frontier_final$cost[i] - frontier_final$cost[i - 1]
      d_qaly <- frontier_final$qaly[i] - frontier_final$qaly[i - 1]
      df$sequential_icer[idx] <- if (d_qaly != 0) d_cost / d_qaly else NA_real_
    }
  }

  return(df)
}

# -----------------------------------------------------------------------------
# Probabilistic sensitivity analysis (PSA)
#
# Runs n_sim Monte Carlo simulations. In each simulation:
#   1. Markov parameters (costs, QALYs, transition probabilities) are drawn
#      from their distributions via sample_params().
#   2. All 43 strategies are run with fixed initial conditions and their
#      costs/QALYs recorded.
#
# Note: initial state distributions are fixed from the Excel decision tree outputs.
# eff_active_* parameters are re-sampled per simulation via beta distributions in config.csv;
# background LTBI prevalence is held fixed (addressed in scenario analyses).
# -----------------------------------------------------------------------------
run_psa <- function(strategies_base, n_sim = 1000, n_cores = NULL) {
  # Samples Markov parameters and eff_active_* detection rates from config.csv distributions.
  # Initial LTBI state distributions are fixed; their uncertainty is addressed in DSA.
  #
  # Parallelisation: uses parallel::mclapply (fork-based; works on macOS/Linux).
  # Each simulation gets a pre-assigned seed derived from the master seed so results
  # are reproducible regardless of core count. n_cores defaults to (detected - 1).

  if (is.null(n_cores)) n_cores <- max(1L, detectCores() - 1L)
  cat(sprintf("\nRunning PSA with %d simulations (%d cores)...\n", n_sim, n_cores))

  # Pre-generate per-simulation seeds from the current RNG state so parallelism
  # does not break reproducibility — same seeds regardless of core count.
  sim_seeds <- sample.int(.Machine$integer.max, n_sim)

  one_sim <- function(i) {
    set.seed(sim_seeds[i])
    sampled <- sample_params()

    # Sample active TB detection efficiency parameters for this PSA simulation.
    # Each eff_active_* parameter has a beta distribution in config.csv.
    eff_sampled <- setNames(
      rbeta(length(eff_active_param_names),
            shape1 = eff_active_dist_p1[eff_active_param_names],
            shape2 = eff_active_dist_p2[eff_active_param_names]),
      eff_active_param_names
    )

    rows    <- vector("list", length(strategies_base))
    for (k in seq_along(strategies_base)) {
      s  <- strategies_base[[k]]
      ek <- s$eff_active_key
      if (!is.null(ek) && !is.na(ek)) {
        # Rescale ActiveDiagnosed/ActiveUndiagnosed using the sampled detection efficiency.
        # Total active TB at entry = base_tp + base_fn (fixed from Excel/prevalence).
        new_eff <- eff_sampled[[ek]]
        s$init["ActiveDiagnosed"]   <- n_c * prev_active * new_eff
        s$init["ActiveUndiagnosed"] <- n_c * prev_active * (1 - new_eff)
      }
      result <- run_strategy(s, params = sampled)
      rows[[k]] <- tibble(
        sim      = i,
        strategy = strategies_base[[k]]$name,
        cost     = result$total_cost,
        qaly     = result$total_qaly
      )
    }
    rows
  }

  raw <- mclapply(seq_len(n_sim), one_sim, mc.cores = n_cores)

  psa_df <- bind_rows(lapply(raw, bind_rows))
  return(psa_df)
}

#------------------------------------------------------------------------------#
# PSA Visualisation: Cost-Effectiveness Plane
#------------------------------------------------------------------------------#
plot_ce_plane <- function(psa_df, reference = "No screening", wtp = 20000) {
  ref_data <- psa_df %>% filter(strategy == reference)

  ce_df <- psa_df %>%
    filter(strategy != reference) %>%
    left_join(ref_data %>% select(sim, ref_cost = cost, ref_qaly = qaly), by = "sim") %>%
    mutate(
      inc_cost = cost - ref_cost,
      inc_qaly = qaly - ref_qaly
    )

  p <- ggplot(ce_df, aes(x = inc_qaly, y = inc_cost, color = strategy)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(intercept = 0, slope = wtp, linetype = "dotted", color = "#7a1a4a", linewidth = 1) +
    annotate("text", x = max(ce_df$inc_qaly) * 0.8, y = wtp * max(ce_df$inc_qaly) * 0.8,
             label = paste0("WTP = £", format(wtp, big.mark = ",")), color = "#7a1a4a", size = 3.5) +
    scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    scale_y_continuous(labels = function(x) paste0("\u00a3", format(round(x), big.mark = ",", scientific = FALSE))) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Incremental QALYs",
      y = "Incremental Cost (\u00a3)",
      title = NULL,
      subtitle = paste0("vs ", reference, " | WTP threshold = £", format(wtp, big.mark = ","), "/QALY"),
      color = "Strategy"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  return(p)
}

#------------------------------------------------------------------------------#
# PSA Visualisation: Cost-Effectiveness Acceptability Curve (CEAC)
#------------------------------------------------------------------------------#
plot_ceac <- function(psa_df, reference = "No screening",
                      wtp_range = seq(0, 50000, by = 1000)) {
  ref_data <- psa_df %>% filter(strategy == reference)
  n_sim <- max(psa_df$sim)

  strategies_to_compare <- unique(psa_df$strategy)

  ceac_data <- expand.grid(wtp = wtp_range, strategy = strategies_to_compare,
                           stringsAsFactors = FALSE)
  ceac_data$prob_ce <- NA

  for (w in wtp_range) {
    # NMB = QALY × WTP − cost; strategy with highest NMB is optimal
    nmb_df <- psa_df %>%
      mutate(nmb = qaly * w - cost) %>%
      group_by(sim) %>%
      mutate(is_best = nmb == max(nmb)) %>%
      ungroup()

    # Probability of being cost-effective at this WTP
    prob_df <- nmb_df %>%
      group_by(strategy) %>%
      summarise(prob_ce = mean(is_best), .groups = "drop")

    for (s in strategies_to_compare) {
      idx <- which(ceac_data$wtp == w & ceac_data$strategy == s)
      prob_val <- prob_df$prob_ce[prob_df$strategy == s]
      if (length(prob_val) > 0) ceac_data$prob_ce[idx] <- prob_val
    }
  }

  # Distinct colours per frontier strategy for CEAC readability.
  ceac_strat_cols <- c(
    "No screening"               = "#2d0a2e",
    "Cough+CXR (TB sx)"          = "#c84b6a",
    "Symptom screen+CXR"         = "#9b5fc0",
    "Parallel Sx+QFT (Ultra)"    = "#5a1a7a"
  )
  # Fall back to auto palette for any strategies not in the named set
  all_ceac_strats <- unique(ceac_data$strategy)
  extra_strats    <- setdiff(all_ceac_strats, names(ceac_strat_cols))
  if (length(extra_strats) > 0) {
    extra_cols <- setNames(project_pal(length(extra_strats)), extra_strats)
    ceac_strat_cols <- c(ceac_strat_cols, extra_cols)
  }

  p <- ggplot(ceac_data, aes(x = wtp, y = prob_ce, color = strategy)) +
    geom_line(linewidth = 1.4) +
    geom_vline(xintercept = 25000, linetype = "dashed", color = "grey40") +
    annotate("text", x = 25000, y = 0.05, label = "NICE\n£25k", color = "grey40", size = 3, hjust = -0.1) +
    geom_vline(xintercept = 35000, linetype = "dashed", color = "grey60") +
    annotate("text", x = 35000, y = 0.05, label = "NICE\n£35k", color = "grey60", size = 3, hjust = -0.1) +
    scale_x_continuous(labels = function(x) paste0("\u00a3", format(x, big.mark = ",", scientific = FALSE))) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_colour_manual(values = ceac_strat_cols) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Willingness-to-Pay Threshold (£/QALY)",
      y = "Probability Cost-Effective",
      title = NULL,
      subtitle = "Probability each strategy is optimal at given WTP thresholds",
      color = "Strategy"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  return(p)
}

# =============================================================================
# BASE CASE ANALYSIS
# Runs the Markov model using central (mean) parameter values from config.csv
# for the no-screening scenario. This produces the reference trajectory
# against which all screening strategies are compared.
# =============================================================================

results <- model(paramsData)
saveRDS(list(state_mat = as.matrix(results$state_membership), n_t = n_t),
        "output/no_screening_state_mat.rds")

# =============================================================================
# BASE CASE RESULTS AND VISUALISATION
# =============================================================================

library(scales)

# -------------------- Summary Output -------------------------------------------
cat("\n")
cat("================================================================================\n")
cat("                    TB SCREENING MODEL - RESULTS SUMMARY                        \n")
cat("================================================================================\n")
cat("\n")
cat("INITIAL STATE DISTRIBUTION:\n")
cat(sprintf("  Uninfected:          %s (%.1f%%)\n",
    format(results$state_membership[1, "Uninfected"], big.mark=","),
    100 * results$state_membership[1, "Uninfected"] / n_c))
cat(sprintf("  Latent Undiagnosed:  %s (%.1f%%)\n",
    format(results$state_membership[1, "LatentUndiagnosed"], big.mark=","),
    100 * results$state_membership[1, "LatentUndiagnosed"] / n_c))
cat(sprintf("  Active Undiagnosed:  %s (%.1f%%)\n",
    format(results$state_membership[1, "ActiveUndiagnosed"], big.mark=","),
    100 * results$state_membership[1, "ActiveUndiagnosed"] / n_c))
cat("\n")
cat(sprintf("OUTCOMES OVER %d YEARS (with half-cycle correction):\n", n_t %/% 12))
cat(sprintf("  Total discounted costs:  £%s\n", format(round(results$total_cost), big.mark=",")))
cat(sprintf("  Total discounted QALYs:  %s\n", format(round(results$total_qaly, 1), big.mark=",")))
cat(sprintf("  Deaths:                  %s\n", format(round(results$state_membership[n_t, "Dead"]), big.mark=",")))
cat("\n")
cat("================================================================================\n")
cat("\n")

# -------------------- Custom Color Palette -------------------------------------
state_colors <- c(
  # Uninfected — light lavender
  "Uninfected"         = "#d8c8e8",
  # Latent TB — purple family (light to dark)
  "LatentUndiagnosed"  = "#b898d0",
  "LatentDiagnosed"    = "#9870b8",
  "LatentTreated"      = "#6040a0",
  "LatentNotreated"    = "#7a50a8",
  "LatentCompleted"    = "#c8b0dc",
  "LatentDiscontinued" = "#8060b0",
  "LatentLtfu"         = "#a888c8",
  # Active TB — rose/pink family (light to dark)
  "ActiveUndiagnosed"  = "#f0a8b8",
  "ActiveDiagnosed"    = "#e06888",
  "ActiveTreated"      = "#800838",
  "ActiveNotreated"    = "#b82858",
  "ActiveCompleted"    = "#f8ccd8",
  "ActiveDiscontinued" = "#982050",
  "ActiveLtfu"         = "#d04870",
  # Dead
  "Dead"               = "#2d2d2d"
)

# -------------------- State Trajectories Plot ----------------------------------
state_mat <- as.matrix(results$state_membership)
state_df <- as.data.frame(state_mat)
state_df$month <- 1:nrow(state_df)
long <- pivot_longer(state_df, cols = -month, names_to = "state", values_to = "n")

state_level_order <- c(
  "Uninfected",
  "LatentUndiagnosed", "LatentDiagnosed", "LatentTreated", "LatentCompleted",
  "LatentNotreated", "LatentDiscontinued", "LatentLtfu",
  "ActiveUndiagnosed", "ActiveDiagnosed", "ActiveTreated", "ActiveCompleted",
  "ActiveNotreated", "ActiveDiscontinued", "ActiveLtfu",
  "Dead"
)
long$state <- factor(long$state, levels = state_level_order)

state_labels <- c(
  "Uninfected"         = "Uninfected",
  "LatentUndiagnosed"  = "LTBI: undiagnosed",
  "LatentDiagnosed"    = "LTBI: diagnosed",
  "LatentTreated"      = "LTBI: on treatment",
  "LatentCompleted"    = "LTBI: treatment completed",
  "LatentNotreated"    = "LTBI: not treated",
  "LatentDiscontinued" = "LTBI: treatment discontinued",
  "LatentLtfu"         = "LTBI: lost to follow-up",
  "ActiveUndiagnosed"  = "Active TB: undiagnosed",
  "ActiveDiagnosed"    = "Active TB: diagnosed",
  "ActiveTreated"      = "Active TB: on treatment",
  "ActiveCompleted"    = "Active TB: treatment completed",
  "ActiveNotreated"    = "Active TB: not treated",
  "ActiveDiscontinued" = "Active TB: treatment discontinued",
  "ActiveLtfu"         = "Active TB: lost to follow-up",
  "Dead"               = "Dead"
)

p_states <- ggplot(long %>% filter(n > 0), aes(x = month, y = n, color = state)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = state_colors,
                     breaks = names(state_labels),
                     labels = state_labels) +
  scale_y_log10(labels = comma, breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  scale_x_continuous(breaks = seq(0, n_t, by = 60), labels = function(x) paste0(x/12, "y")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Time (years)",
    y = "Number of individuals (log\u2081\u2080 scale)",
    color = "Health state"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10)
  )

# -------------------- Display & Save State Trajectories -----------------------
print(p_states)
ggsave("output/state_trajectories.png", p_states, width = 12, height = 7, dpi = 300)
cat("Plots saved to output/ folder\n")

# ================================================================================
# EPIDEMIOLOGY & IN-DEPTH ANALYSIS PLOTS
# ================================================================================

# -------------------- 1. TB Incidence Plot -------------------------------------
# Calculates monthly active TB incidence by counting new transitions from
# latent states into ActiveUndiagnosed. This reflects genuine disease
# progression from latent to active TB rather than the snapshot state counts.
# Re-entry from ActiveDiscontinued/ActiveLtfu is excluded as these individuals
# already had active TB and are not new incident cases.
# Sources of new active cases per cycle:
#   LatentUndiagnosed -> ActiveUndiagnosed (observed incidence; Berrocal-Almanza 2022)
#   LatentDiagnosed   -> ActiveUndiagnosed (progression despite diagnosis)
#   LatentNotreated   -> ActiveUndiagnosed (untreated LTBI progression)
#   LatentDiscontinued -> ActiveUndiagnosed (partial treatment failure)
#   LatentLtfu        -> ActiveUndiagnosed (lost to follow-up; no treatment protection)
active_states <- c("ActiveUndiagnosed", "ActiveDiagnosed", "ActiveTreated",
                   "ActiveNotreated", "ActiveCompleted", "ActiveDiscontinued", "ActiveLtfu")

m_p_base <- results$m_p
new_active_cases <- numeric(n_t)
new_active_cases[1] <- 0
for (t in 2:n_t) {
  # People transitioning INTO ActiveUndiagnosed from latent/other states
  new_active_cases[t] <-
    state_mat[t-1, "LatentUndiagnosed"] * m_p_base["LatentUndiagnosed", "ActiveUndiagnosed"] +
    state_mat[t-1, "LatentDiagnosed"] * m_p_base["LatentDiagnosed", "ActiveUndiagnosed"] +
    state_mat[t-1, "LatentNotreated"] * m_p_base["LatentNotreated", "ActiveUndiagnosed"] +
    state_mat[t-1, "LatentDiscontinued"] * m_p_base["LatentDiscontinued", "ActiveUndiagnosed"] +
    state_mat[t-1, "LatentLtfu"] * m_p_base["LatentLtfu", "ActiveUndiagnosed"]
}

person_months <- rowSums(state_mat[, v_state_names[v_state_names != "Dead"]])
incidence_per_1000_py <- (new_active_cases / person_months) * 1000 * 12

df_incidence <- tibble(month = 1:n_t, incidence = incidence_per_1000_py)

p_incidence <- ggplot(df_incidence, aes(x = month, y = incidence)) +
  geom_line(color = "#c84b6a", linewidth = 1.2) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2, color = "#c84b6a", fill = "#f0c0cc") +
  scale_x_continuous(breaks = seq(0, n_t, by = 24), labels = function(x) paste0(x/12, "y")) +
  theme_minimal(base_size = 14) +
  labs(x = "Time", y = "Incidence per 1,000 person-years",
       title = NULL,
       subtitle = "No Screening scenario | Raw per-cycle rate (line) + LOESS smooth (band)") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank())

# -------------------- 2. TB Prevalence Plot ------------------------------------
latent_states <- c("LatentUndiagnosed", "LatentDiagnosed", "LatentTreated",
                   "LatentNotreated", "LatentCompleted", "LatentDiscontinued", "LatentLtfu")

df_prevalence <- tibble(
  month = rep(1:n_t, 3),
  category = factor(rep(c("Uninfected", "Latent TB", "Active TB"), each = n_t),
                    levels = c("Active TB", "Latent TB", "Uninfected")),
  count = c(state_mat[, "Uninfected"], rowSums(state_mat[, latent_states]), rowSums(state_mat[, active_states]))
)

p_prevalence <- ggplot(df_prevalence, aes(x = month, y = count, fill = category)) +
  geom_area(alpha = 0.8) +
  scale_fill_manual(values = c("Active TB" = "#c84b6a", "Latent TB" = "#9870b8", "Uninfected" = "#f0e8f8")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(0, n_t, by = 24), labels = function(x) paste0(x/12, "y")) +
  theme_minimal(base_size = 14) +
  labs(x = "Time", y = "Number of individuals",
       title = NULL,
       subtitle = "Distribution of uninfected, latent, and active TB", fill = "Status") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(), legend.position = "right")

# -------------------- 3. Treatment Cascade (Cumulative Flows) ------------------
# The treatment cascade quantifies how many individuals pass through each
# stage of the care pathway over 20 years. Cumulative flows are estimated
# by summing the transition probabilities applied to state membership at
# each cycle — this gives true individual counts rather than person-months.
# Initial state values (month 1) are added to capture those detected at entry.
m_p_base <- results$m_p

# For treatment cascade and cumulative outcomes: use Parallel Sx+QFT (Ultra)
# — the frontier parallel IGRA strategy (seq. ICER £7,745/QALY). This shows
# the full LTBI treatment cascade (11,526 LTBI detected via specificity-corrected QFT-GIT).
# Inline init: TP=874.6, FN=125.4 per 100k (decision tree; Parallel Sx+QFT Ultra);
# LTBI = n_ltbi × igra_sens_qft − FP (specificity-corrected; Pai 2008)
init_parallel_sx_qft <- setNames(rep(0, length(v_state_names)), v_state_names)
ltbi_det_pcqft <- max(0, n_c * prev_ltbi * igra_sens_qft - (1 - igra_spec_qft) * n_c * prev_uninfected)  # specificity-corrected: ~11,526
init_parallel_sx_qft["ActiveDiagnosed"]   <- 874.6 * (n_c / excel_denom)
init_parallel_sx_qft["ActiveUndiagnosed"] <- 125.4 * (n_c / excel_denom)
init_parallel_sx_qft["LatentDiagnosed"]   <- ltbi_det_pcqft
init_parallel_sx_qft["LatentUndiagnosed"] <- n_c * prev_ltbi - ltbi_det_pcqft
init_parallel_sx_qft["Uninfected"]        <- n_c * prev_uninfected
results_psqft   <- model(paramsData, init_dist = init_parallel_sx_qft)
state_mat_psqft <- as.matrix(results_psqft$state_membership)
state_mat       <- state_mat_psqft   # override: cascade + cumulative use Parallel Sx+QFT (Ultra)

# Note: the person-month sums below are computed but superseded by the
# transition-flow approach (cum_flow_* variables). They are retained here
# for reference only.
cum_ltbi_diagnosed <- sum(state_mat[, "LatentDiagnosed"]) +
  sum(state_mat[, "LatentTreated"]) + sum(state_mat[, "LatentNotreated"]) +
  sum(state_mat[, "LatentCompleted"]) + sum(state_mat[, "LatentDiscontinued"]) +
  sum(state_mat[, "LatentLtfu"])
cum_ltbi_treated <- sum(state_mat[, "LatentTreated"]) +
  sum(state_mat[, "LatentCompleted"]) + sum(state_mat[, "LatentDiscontinued"]) +
  sum(state_mat[, "LatentLtfu"])
cum_ltbi_completed    <- sum(state_mat[, "LatentCompleted"])
cum_active_diagnosed  <- sum(state_mat[, "ActiveDiagnosed"]) +
  sum(state_mat[, "ActiveTreated"]) + sum(state_mat[, "ActiveNotreated"]) +
  sum(state_mat[, "ActiveCompleted"]) + sum(state_mat[, "ActiveDiscontinued"]) +
  sum(state_mat[, "ActiveLtfu"])
cum_active_treated    <- sum(state_mat[, "ActiveTreated"]) +
  sum(state_mat[, "ActiveCompleted"]) + sum(state_mat[, "ActiveDiscontinued"]) +
  sum(state_mat[, "ActiveLtfu"])
cum_active_completed  <- sum(state_mat[, "ActiveCompleted"])

# Transition-flow method: count individuals moving between cascade stages
# each month and accumulate over the full time horizon.
cum_flow_ltbi_diag <- 0
cum_flow_ltbi_treat <- 0
cum_flow_ltbi_comp <- 0
cum_flow_active_diag <- 0
cum_flow_active_treat <- 0
cum_flow_active_comp <- 0
for (t in 2:n_t) {
  cum_flow_ltbi_diag  <- cum_flow_ltbi_diag +
    state_mat[t-1, "LatentUndiagnosed"] * m_p_base["LatentUndiagnosed", "LatentDiagnosed"]
  cum_flow_ltbi_treat <- cum_flow_ltbi_treat +
    state_mat[t-1, "LatentDiagnosed"] * m_p_base["LatentDiagnosed", "LatentTreated"]
  cum_flow_ltbi_comp  <- cum_flow_ltbi_comp +
    state_mat[t-1, "LatentTreated"] * m_p_base["LatentTreated", "LatentCompleted"]
  cum_flow_active_diag  <- cum_flow_active_diag +
    state_mat[t-1, "ActiveUndiagnosed"] * m_p_base["ActiveUndiagnosed", "ActiveDiagnosed"]
  cum_flow_active_treat <- cum_flow_active_treat +
    state_mat[t-1, "ActiveDiagnosed"] * m_p_base["ActiveDiagnosed", "ActiveTreated"]
  cum_flow_active_comp  <- cum_flow_active_comp +
    state_mat[t-1, "ActiveTreated"] * m_p_base["ActiveTreated", "ActiveCompleted"]
}
# Add initial state counts (people who start already diagnosed/treated)
cum_flow_ltbi_diag  <- cum_flow_ltbi_diag + state_mat[1, "LatentDiagnosed"]
cum_flow_active_diag <- cum_flow_active_diag + state_mat[1, "ActiveDiagnosed"]

cascade_data <- tibble(
  stage = factor(c("LTBI Diagnosed", "LTBI Treated", "LTBI Completed",
                   "Active Diagnosed", "Active Treated", "Active Completed"),
                 levels = c("LTBI Diagnosed", "LTBI Treated", "LTBI Completed",
                           "Active Diagnosed", "Active Treated", "Active Completed")),
  count = c(cum_flow_ltbi_diag, cum_flow_ltbi_treat, cum_flow_ltbi_comp,
            cum_flow_active_diag, cum_flow_active_treat, cum_flow_active_comp),
  type = factor(c(rep("Latent TB", 3), rep("Active TB", 3)))
)

p_cascade <- ggplot(cascade_data, aes(x = stage, y = count, fill = type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = round(count)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Latent TB" = "#9870b8", "Active TB" = "#c84b6a")) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "Number of individuals",
       title = NULL,
       subtitle = sprintf("Parallel Sx+QFT (Ultra) strategy | Cumulative individuals passing through each care stage over %d years", n_t %/% 12), fill = "") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

# -------------------- 4. Cumulative Outcomes -----------------------------------
df_cumulative <- tibble(
  month = rep(1:n_t, 4),
  outcome = factor(rep(c("Deaths", "Treatment Completed", "Lost to Follow-up", "Treatment Discontinued"), each = n_t),
                   levels = c("Deaths", "Treatment Completed", "Lost to Follow-up", "Treatment Discontinued")),
  count = c(state_mat[, "Dead"],
            state_mat[, "LatentCompleted"] + state_mat[, "ActiveCompleted"],
            state_mat[, "LatentLtfu"] + state_mat[, "ActiveLtfu"],
            state_mat[, "LatentDiscontinued"] + state_mat[, "ActiveDiscontinued"])
)

p_cumulative <- ggplot(df_cumulative, aes(x = month, y = count, color = outcome)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Deaths" = "#2d0a2e", "Treatment Completed" = "#c84b6a",
                                 "Lost to Follow-up" = "#9870b8", "Treatment Discontinued" = "#e8a0b8")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(0, n_t, by = 24), labels = function(x) paste0(x/12, "y")) +
  theme_minimal(base_size = 14) +
  labs(x = "Time", y = "Cumulative count",
       title = NULL,
       subtitle = "Parallel Sx+QFT (Ultra) strategy | Deaths, treatment completions, LTFU, and discontinuations", color = "Outcome") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(), legend.position = "right")

# -------------------- 5. Stacked State Distribution ----------------------------
state_categories <- tibble(
  state = v_state_names,
  category = case_when(
    state == "Uninfected" ~ "Uninfected",
    state == "Dead" ~ "Dead",
    grepl("^Latent", state) ~ "Latent TB",
    grepl("^Active", state) ~ "Active TB"
  )
)

saveRDS(list(state_mat = state_mat, n_t = n_t, n_c = n_c,
             state_categories = state_categories),
        "output/pie_cache.rds")

df_stacked <- as.data.frame(state_mat) %>%
  mutate(month = 1:n_t) %>%
  pivot_longer(cols = -month, names_to = "state", values_to = "n") %>%
  left_join(state_categories, by = "state") %>%
  group_by(month, category) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(category = factor(category, levels = c("Dead", "Active TB", "Latent TB", "Uninfected")))

p_stacked <- ggplot(df_stacked, aes(x = month, y = n, fill = category)) +
  geom_area(alpha = 0.9) +
  scale_fill_manual(values = c("Dead" = "#2d2d2d", "Active TB" = "#c84b6a",
                                "Latent TB" = "#9870b8", "Uninfected" = "#f0e8f8")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(0, n_t, by = 24), labels = function(x) paste0(x/12, "y")) +
  theme_minimal(base_size = 14) +
  labs(x = "Time", y = "Number of individuals",
       title = NULL,
       subtitle = "Parallel Sx+QFT (Ultra) strategy | Stacked area chart showing disease progression", fill = "Category") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(), legend.position = "right")

# -------------------- 6. Final State Pie Chart ---------------------------------
final_dist <- df_stacked %>%
  filter(month == n_t) %>%
  mutate(pct = n / sum(n) * 100, label = paste0(category, "\n", round(pct, 1), "%"))

final_dist <- final_dist %>%
  mutate(n_fmt = format(round(n), big.mark = ","))

# Pie chart: ggplot2 coord_polar stacks in REVERSE factor level order, so
# midpoints must be computed over data sorted in reverse factor level order.
pie_cols <- c("Dead" = "#7a1a4a", "Active TB" = "#c84b6a",
              "Latent TB" = "#e8a0b8", "Uninfected" = "#fce8f0")
final_dist_pie <- final_dist %>%
  mutate(category = factor(category,
           levels = c("Dead", "Active TB", "Latent TB", "Uninfected"))) %>%
  arrange(desc(category)) %>%           # reverse order for correct cumsum
  mutate(midpoint   = cumsum(n) - n / 2,
         label_text = paste0(category, "\n", n_fmt, "\n(", round(pct, 1), "%)"))
p_pie <- ggplot(final_dist_pie, aes(x = "", y = n, fill = category)) +
  geom_col(width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = pie_cols) +
  geom_label(aes(y = midpoint, label = label_text),
             x = 1.6, size = 2.8, show.legend = FALSE,
             label.padding = unit(0.15, "lines"), linewidth = 0.3) +
  theme_void(base_size = 14) +
  labs(title = paste0("Final Population Distribution (Year ", n_t %/% 12, ")"),
       subtitle = paste0("Parallel Sx+QFT (Ultra) strategy | Total population: ", prettyNum(n_c, big.mark = ",")),
       fill = "Category") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(color = "grey40", hjust = 0.5), legend.position = "right")

# -------------------- 7. Cost Breakdown Bar Chart ------------------------------
cost_by_state <- tibble(
  state = v_state_names,
  total_cost = colSums(state_mat * rep(as.numeric(results$m_payoffsCost), each = n_t))
) %>%
  left_join(state_categories, by = "state") %>%
  group_by(category) %>%
  summarise(total_cost = sum(total_cost), .groups = "drop") %>%
  mutate(category = factor(category, levels = c("Active TB", "Latent TB", "Uninfected", "Dead")),
         pct = total_cost / sum(total_cost) * 100)

p_cost_breakdown <- ggplot(cost_by_state, aes(x = reorder(category, -total_cost), y = total_cost, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("£", format(round(total_cost), big.mark = ","), "\n(", round(pct, 1), "%)")),
            vjust = -0.2, size = 4) +
  scale_fill_manual(values = c("Dead" = "#2d2d2d", "Active TB" = "#c84b6a",
                                "Latent TB" = "#9870b8", "Uninfected" = "#f0e8f8")) +
  scale_y_continuous(labels = label_comma(prefix = "£"), expand = expansion(mult = c(0, 0.2))) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "Total undiscounted costs",
       title = NULL,
       subtitle = "Parallel Sx+QFT (Ultra) strategy | Undiscounted | Which disease states drive costs?") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(), legend.position = "none")

# -------------------- 8. Summary Statistics ------------------------------------
cat("\n")
cat("================================================================================\n")
cat("                    EPIDEMIOLOGY SUMMARY STATISTICS                             \n")
cat("================================================================================\n")
cat("\n")
cat("DISEASE BURDEN:\n")
cat(sprintf("  Peak LTBI prevalence:     %s (%.1f%% of cohort)\n",
    format(round(max(rowSums(state_mat[, latent_states]))), big.mark = ","),
    100 * max(rowSums(state_mat[, latent_states])) / n_c))
cat(sprintf("  Peak Active TB:           %s (%.2f%% of cohort)\n",
    format(round(max(rowSums(state_mat[, active_states]))), big.mark = ","),
    100 * max(rowSums(state_mat[, active_states])) / n_c))
cat(sprintf("  Total deaths:             %s (%.2f%% mortality)\n",
    format(round(state_mat[n_t, "Dead"]), big.mark = ","),
    100 * state_mat[n_t, "Dead"] / n_c))
cat("\nTREATMENT OUTCOMES:\n")
cat(sprintf("  LTBI treatment completed: %s\n", format(round(state_mat[n_t, "LatentCompleted"]), big.mark = ",")))
cat(sprintf("  Active TB cured:          %s\n", format(round(state_mat[n_t, "ActiveCompleted"]), big.mark = ",")))
cat(sprintf("  Lost to follow-up:        %s\n",
    format(round(state_mat[n_t, "LatentLtfu"] + state_mat[n_t, "ActiveLtfu"]), big.mark = ",")))
cat("\nECONOMIC ANALYSIS (with half-cycle correction):\n")
total_costs <- results$total_cost
total_qalys <- results$total_qaly
cat(sprintf("  Total discounted costs:   £%s\n", format(round(total_costs), big.mark = ",")))
cat(sprintf("  Total discounted QALYs:   %s\n", format(round(total_qalys, 1), big.mark = ",")))
cat(sprintf("  Cost per person:          £%s\n", format(round(total_costs / n_c), big.mark = ",")))
cat(sprintf("  QALYs per person:         %.2f\n", total_qalys / n_c))
cat("\n================================================================================\n\n")

# -------------------- Save Epidemiology Plots ----------------------------------
ggsave("output/tb_incidence.png", p_incidence, width = 10, height = 5, dpi = 300)
ggsave("output/tb_prevalence.png", p_prevalence, width = 10, height = 6, dpi = 300)
ggsave("output/treatment_cascade.png", p_cascade, width = 10, height = 6, dpi = 300)
ggsave("output/cumulative_outcomes.png", p_cumulative, width = 10, height = 6, dpi = 300)
ggsave("output/state_distribution_stacked.png", p_stacked, width = 10, height = 6, dpi = 300)
ggsave("output/final_state_pie.png", p_pie, width = 8, height = 6, dpi = 300)
ggsave("output/cost_breakdown.png", p_cost_breakdown, width = 8, height = 6, dpi = 300)

cat("Epidemiology plots saved to output/ folder\n")

# =============================================================================
# DIAGNOSTIC STRATEGY COMPARISON — BASE CASE ICER TABLE
# Runs all 43 strategies using central parameter values.
# Results are compared against no screening and ranked by ICER.
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("                    DIAGNOSTIC STRATEGY COMPARISON                              \n")
cat("================================================================================\n\n")

strategies <- create_diagnostic_strategies()

# Run all 43 strategies with base case (mean) parameter values
strategy_results <- mclapply(strategies, run_strategy,
                             mc.cores = max(1L, detectCores() - 1L))

# Calculate ICER table
icer_table <- calculate_icer(strategy_results)

cat("ICER TABLE (sequential incremental analysis with extended dominance):\n")
cat(sprintf("%-22s %11s %10s %11s %11s %14s %15s  %s\n",
            "Strategy", "Cost/person", "QALYs/pp", "Inc.Cost", "Inc.QALY",
            "ICER(vs ref)", "Sequential ICER", "Dominance"))
cat(strrep("-", 110), "\n")
for (i in 1:nrow(icer_table)) {
  icer_str     <- if (is.na(icer_table$icer[i])) "Ref"
                  else paste0("£", format(round(icer_table$icer[i]), big.mark = ","))
  seq_icer_str <- if (is.na(icer_table$sequential_icer[i])) "-"
                  else paste0("£", format(round(icer_table$sequential_icer[i]), big.mark = ","))
  cat(sprintf("%-22s %10s %10.4f %10s %11.4f %14s %15s  %s\n",
              icer_table$strategy[i],
              paste0("£", format(round(icer_table$cost_per_person[i]), big.mark = ",")),
              icer_table$qaly_per_person[i],
              paste0("£", format(round(icer_table$inc_cost[i] / n_c), big.mark = ",")),
              icer_table$inc_qaly[i] / n_c,
              icer_str,
              seq_icer_str,
              icer_table$dominance[i]))
}
cat("\n================================================================================\n")

# -------------------- Efficient frontier summary --------------------------------
frontier_rows <- icer_table[icer_table$dominance %in% c("ref", "non-dominated"), ]
frontier_rows <- frontier_rows[order(frontier_rows$cost), ]

dominated_strategies <- icer_table$strategy[icer_table$dominance == "simply dominated"]
ext_dom_strategies   <- icer_table$strategy[icer_table$dominance == "extendedly dominated"]

cat("\n")
cat("================================================================================\n")
cat("                         EFFICIENT FRONTIER SUMMARY                             \n")
cat("================================================================================\n\n")
cat("COST-EFFECTIVE FRONTIER (strategies on the non-dominated frontier):\n")
cat(sprintf("  %-24s  Cost/pp  Seq. ICER vs previous\n", "Strategy"))
cat(strrep("-", 65), "\n")
for (i in seq_len(nrow(frontier_rows))) {
  seq_str <- if (is.na(frontier_rows$sequential_icer[i])) "Reference"
             else paste0("£", format(round(frontier_rows$sequential_icer[i]), big.mark = ","), "/QALY")
  cat(sprintf("  %-24s  £%-7s  %s\n",
              frontier_rows$strategy[i],
              format(round(frontier_rows$cost_per_person[i]), big.mark = ","),
              seq_str))
}
cat("\n")

if (length(dominated_strategies) > 0) {
  cat(sprintf("SIMPLY DOMINATED (%d strategies — higher cost AND lower QALYs than a cheaper alternative):\n",
              length(dominated_strategies)))
  for (s in dominated_strategies) cat(sprintf("  - %s\n", s))
  cat("\n")
}
if (length(ext_dom_strategies) > 0) {
  cat(sprintf("EXTENDEDLY DOMINATED (%d strategies — worse value than linear interpolation on frontier):\n",
              length(ext_dom_strategies)))
  for (s in ext_dom_strategies) cat(sprintf("  - %s\n", s))
  cat("\n")
}

cat("INTERPRETATION (NICE WTP threshold: £25,000/QALY; updated Dec 2025, effective Apr 2026):\n")
for (i in seq_len(nrow(frontier_rows))) {
  if (is.na(frontier_rows$sequential_icer[i])) next  # skip reference
  wtp <- 25000
  ce_status <- if (frontier_rows$sequential_icer[i] <= wtp) "COST-EFFECTIVE" else "NOT cost-effective at £25k"
  cat(sprintf("  %-24s  seq. ICER £%s → %s\n",
              frontier_rows$strategy[i],
              format(round(frontier_rows$sequential_icer[i]), big.mark = ","),
              ce_status))
}
cat("\n================================================================================\n\n")

# =============================================================================
# EVENT COUNTS — FRONTIER STRATEGIES VS NO SCREENING
#
# Reports the key mechanism statistics for each frontier strategy:
#   - LTBI detected at entry (from init vector — the structural gap between
#     sequential strategies (~100/100k) and parallel IGRA (~11,526 (QFT) / 9,980 (T-SPOT) per 100k))
#   - Active TB detected at entry (TP from Zenner 2025 decision tree)
#   - Deaths over 55-year lifetime horizon (undiscounted; from Markov Dead state)
#   - Deaths averted vs No Screening
#
# =============================================================================
no_screen_idx    <- which(sapply(strategy_results, function(x) x$strategy_name == "No screening"))
no_screen_res    <- strategy_results[[no_screen_idx]]
no_screen_deaths <- no_screen_res$state_membership[nrow(no_screen_res$state_membership), "Dead"]
no_screen_ltbi   <- unname(strategies[["No_screening"]]$init["LatentDiagnosed"])     # 0
no_screen_atb    <- unname(strategies[["No_screening"]]$init["ActiveDiagnosed"])     # 0

cat("================================================================================\n")
cat("               EVENT COUNTS — FRONTIER STRATEGIES VS NO SCREENING               \n")
cat("================================================================================\n\n")
cat(sprintf("  %-34s  %10s  %10s  %10s  %12s\n",
            "Strategy", "LTBI det.", "ATB det.", "Deaths", "Deaths avrt"))
cat(strrep("-", 82), "\n")

for (i in seq_len(nrow(frontier_rows))) {
  strat_nm  <- frontier_rows$strategy[i]
  strat_key <- names(strategies)[sapply(names(strategies), function(k) strategies[[k]]$name == strat_nm)]
  if (length(strat_key) == 0) next
  sr_idx    <- which(sapply(strategy_results, function(x) x$strategy_name == strat_nm))
  if (length(sr_idx) == 0) next
  sr        <- strategy_results[[sr_idx]]
  deaths_i  <- sr$state_membership[nrow(sr$state_membership), "Dead"]
  ltbi_det  <- unname(strategies[[strat_key[1]]]$init["LatentDiagnosed"])
  atb_det   <- unname(strategies[[strat_key[1]]]$init["ActiveDiagnosed"])
  cat(sprintf("  %-34s  %10.0f  %10.0f  %10.0f  %12.0f\n",
              strat_nm,
              ltbi_det, atb_det, deaths_i,
              no_screen_deaths - deaths_i))
}
cat("\n")
cat("  Notes: LTBI detected = persons in LatentDiagnosed state at t=0 (entry screening)\n")
cat("         ATB detected  = true positives (TP) from Zenner 2025 decision tree\n")
cat("         Deaths        = undiscounted cumulative deaths over 55-year horizon\n\n")

# Save event counts for all strategies
event_counts_df <- bind_rows(lapply(names(strategies), function(k) {
  sr_idx <- which(sapply(strategy_results, function(x) x$strategy_name == strategies[[k]]$name))
  if (length(sr_idx) == 0) return(NULL)
  sr <- strategy_results[[sr_idx]]
  tibble(
    strategy        = strategies[[k]]$name,
    ltbi_detected   = unname(strategies[[k]]$init["LatentDiagnosed"]),
    active_detected = unname(strategies[[k]]$init["ActiveDiagnosed"]),
    deaths_55yr     = sr$state_membership[nrow(sr$state_membership), "Dead"],
    deaths_averted  = no_screen_deaths - sr$state_membership[nrow(sr$state_membership), "Dead"]
  )
}))
write.csv(event_counts_df, "output/csv/event_counts_basecase.csv", row.names = FALSE)
cat("Saved: output/csv/event_counts_basecase.csv\n\n")

# Save base case ICER table as CSV (all 43 strategies), with detection columns added
write.csv(icer_table %>%
  mutate(cost_per_person = cost / n_c, qaly_per_person = qaly / n_c) %>%
  left_join(event_counts_df %>% select(strategy, ltbi_detected, active_detected),
            by = "strategy"),
  "output/csv/icer_table_basecase.csv", row.names = FALSE)
cat("Saved: output/csv/icer_table_basecase.csv\n\n")

# Non-dominated strategy names (ref + efficient frontier) — used to focus labels
# and filter PSA plots. 43 strategies make fully-labelled/coloured plots unreadable.
non_dominated_names <- icer_table$strategy[
  icer_table$dominance %in% c("ref", "non-dominated")
]

# -------------------- Base case cost-effectiveness scatter plot ----------------
# Each point is coloured by dominance status:
#   dark rose = reference (No screening) or efficient frontier
#   light pink = dominated (simply or extendedly)
# No WTP line — this is total cost vs total QALYs (not incremental); WTP line
# appears in the incremental CE plane (plot_ce_plane) and ICER forest plot.
icer_table$ce_status <- dplyr::case_when(
  icer_table$dominance == "ref"                  ~ "Reference",
  icer_table$dominance == "non-dominated"        ~ "Efficient frontier",
  icer_table$dominance == "extendedly dominated" ~ "Extendedly dominated",
  TRUE                                           ~ "Dominated"
)
# Shape aesthetic: triangle for parallel IGRA strategies, circle for all others
icer_table$strategy_type <- ifelse(
  grepl("^Parallel", icer_table$strategy, ignore.case = FALSE),
  "Parallel IGRA", "Sequential"
)

p_icer <- ggplot(icer_table, aes(x = qaly / n_c, y = cost / n_c,
                                   colour = ce_status, shape = strategy_type)) +
  geom_point(aes(size = strategy_type)) +
  ggrepel::geom_text_repel(
    data   = icer_table %>% filter(dominance %in% c("ref", "non-dominated")),
    aes(label = strategy,
        fontface = ifelse(strategy_type == "Parallel IGRA", "bold", "plain")),
    size = 3, max.overlaps = Inf, colour = "grey20",
    box.padding = 0.8, point.padding = 0.5, show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("Reference" = "#7a1a4a", "Efficient frontier" = "#c84b6a",
               "Extendedly dominated" = "#a0506e", "Dominated" = "#e0c0cc"),
    breaks = c("Reference", "Efficient frontier", "Extendedly dominated", "Dominated"),
    name   = NULL
  ) +
  scale_shape_manual(
    values = c("Parallel IGRA" = 17L, "Sequential" = 16L),
    name   = NULL
  ) +
  scale_size_manual(
    values = c("Parallel IGRA" = 4.5, "Sequential" = 3.5),
    guide  = "none"
  ) +
  theme_minimal(base_size = 14) +
  labs(x = "QALYs per person", y = "Cost per person (\u00a3)",
       title = NULL,
       subtitle = NULL) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave("output/ce_scatter_basecase.png", p_icer, width = 14, height = 8, dpi = 300)

# =============================================================================
# STRATEGY SPACE VISUALISATION
#
# Plots all strategies from Zenner et al. 2025 Eur Respir J decision tree to contextualise the
# strategies selected for the Markov analysis.  Uses excel_ic already loaded
# at the top of the script.
#
# Plot 1 — scatter: TP detected (x) vs diagnostic cost per person (y),
#           coloured by test family; selected strategies labelled.
# Plot 2 — horizontal bar: TP detected, ranked descending.
# =============================================================================
cat("\n")
cat("================================================================================\n")
cat("         STRATEGY SPACE VISUALISATION (all Excel strategies)                   \n")
cat("================================================================================\n\n")

# Build data frame from all rows in excel_ic ----------------------------------------
strategy_space <- data.frame(
  name      = as.character(excel_ic[[1]]),
  TP        = suppressWarnings(as.numeric(excel_ic[[4]])),
  Tot_costs = suppressWarnings(as.numeric(excel_ic[[6]])),
  stringsAsFactors = FALSE
)
strategy_space <- strategy_space[
  !is.na(strategy_space$TP) & !is.na(strategy_space$Tot_costs) &
  strategy_space$name != "" & !is.na(strategy_space$name), ]
strategy_space$cost_per_person <- strategy_space$Tot_costs / 100000

strategy_space$family <- sapply(strategy_space$name, assign_family)

# All feasible strategies (non-NA costs) are included in the Markov analysis
strategy_space$selected <- !is.na(strategy_space$Tot_costs)

strategy_space$display_name <- sapply(strategy_space$name, clean_excel_name)

# --- Plot 1: Scatter — TP detected vs diagnostic cost per person -----------------
# Points coloured by primary screening test family.
family_pal  <- c(
  "CXR / Symptom screen"    = "#e8a0b8",
  "QFT-GIT"                 = "#9b5fc0",
  "T-SPOT.TB"               = "#5a1a7a",
  "TST Mantoux"             = "#a03060",
  "Molecular (Xpert/Ultra)" = "#7a1a4a"
)

p_space_scatter <- ggplot(
  strategy_space,
  aes(x = TP, y = cost_per_person, colour = family)
) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = display_name), size = 2.5,
                  max.overlaps = 20, seed = 42, colour = "#333333",
                  show.legend = FALSE) +
  scale_colour_manual(values = family_pal, name = "Primary screening test") +
  theme_minimal(base_size = 13) +
  labs(
    x     = "True positives detected (active TB, per 100,000 screened)",
    y     = "Diagnostic cost per person (\u00a3)",
    title = NULL
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 15),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

ggsave("output/strategy_space_scatter.png", p_space_scatter,
       width = 12, height = 7, dpi = 300)
cat("Saved: output/strategy_space_scatter.png\n")

# --- Plot 2: Two-panel bar — Active TB detection (left) + LTBI detection (right)
# Uses event_counts_df (all 43 strategies); colour-coded by pathway type.
# Parallel IGRA strategies: medium purple (QFT) / dark purple (T-SPOT) — consistent with CE scatter.
# TST is NOT an IGRA — only QFT-GIT and T-SPOT.TB are IGRAs.
# TST strategies grouped with non-IGRA (CXR/symptom screen).
classify_igra <- function(nm) {
  if (grepl("^Parallel", nm, ignore.case = TRUE)) return("Parallel IGRA (QFT-GIT / T-SPOT.TB)")
  if (grepl("QFT-GIT|T-SPOT", nm, ignore.case = TRUE)) return("IGRA-inclusive (sequential pathway)")
  return("Non-IGRA (CXR / symptom screen / TST)")
}

igra_pal_bar <- c(
  "Non-IGRA (CXR / symptom screen / TST)"  = "#f5d0dc",
  "IGRA-inclusive (sequential pathway)"     = "#c090d0",
  "Parallel IGRA (QFT-GIT / T-SPOT.TB)"    = "#7b2d8b"
)

combined_bar_df <- event_counts_df %>%
  filter(strategy != "No screening") %>%
  mutate(
    igra_type = factor(sapply(strategy, classify_igra),
      levels = c("Non-IGRA (CXR / symptom screen / TST)",
                 "IGRA-inclusive (sequential pathway)",
                 "Parallel IGRA (QFT-GIT / T-SPOT.TB)"))
  ) %>%
  arrange(igra_type, active_detected) %>%
  mutate(strategy_fct = factor(strategy, levels = strategy))

shared_bar_fill <- scale_fill_manual(values = igra_pal_bar, name = NULL)

p_active_bar <- ggplot(combined_bar_df,
    aes(x = active_detected / n_c * 100, y = strategy_fct, fill = igra_type)) +
  geom_col(width = 0.7, colour = NA) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.08)),
    labels = function(x) paste0(x, "%")
  ) +
  shared_bar_fill +
  theme_minimal(base_size = 9.5) +
  labs(
    x     = "Active TB cases detected (% of cohort)",
    y     = NULL,
    title = NULL
  ) +
  theme(
    plot.title         = element_text(face = "bold", size = 11),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 7)
  )

p_ltbi_bar <- ggplot(combined_bar_df,
    aes(x = ltbi_detected / n_c * 100, y = strategy_fct, fill = igra_type)) +
  geom_col(width = 0.7, colour = NA) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    labels = function(x) paste0(x, "%")
  ) +
  shared_bar_fill +
  theme_minimal(base_size = 9.5) +
  labs(
    x     = "LTBI cases detected (% of cohort)",
    y     = NULL,
    title = NULL
  ) +
  theme(
    plot.title         = element_text(face = "bold", size = 11),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank()
  )

p_space_bar <- (p_active_bar + p_ltbi_bar) +
  plot_layout(guides = "collect", widths = c(1.15, 1)) &
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 9),
        legend.key.size  = unit(0.8, "lines"))

p_space_bar <- p_space_bar +
  plot_annotation(
    subtitle = "42 active screening strategies shown (No screening reference excluded) \u00b7 grouped by test type \u00b7 sorted by active TB detection within each group \u00b7 percentage of cohort (n = 100,000)",
    theme = theme(
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

ggsave("output/strategy_space_bar.png", p_space_bar,
       width = 13, height = 12, dpi = 300)
cat("Saved: output/strategy_space_bar.png\n")

# =============================================================================
# LTBI DETECTION GAP FIGURE
#
# Horizontal bar chart showing LTBI detected at entry per 100,000 across all
# 43 strategies, coloured by pathway type. This is the central structural
# finding: sequential strategies detect ~72–116 LTBI/100k (driven by the
# decision tree's sequential pathway design — IGRA applied only to
# symptom-screen positives); parallel IGRA strategies detect 14,774–15,664/100k
# (IGRA applied to all migrants). The ~100× gap is the mechanism behind the
# frontier composition.
# =============================================================================

# Classify each strategy by pathway type for colouring
classify_pathway <- function(nm) {
  if (nm == "No screening")                          return("No screening")
  if (grepl("^Parallel.*QFT",  nm, ignore.case = TRUE)) return("Parallel IGRA (QFT)")
  if (grepl("^Parallel.*T-SPOT|^Parallel.*TSPT|^Parallel.*tspt", nm, ignore.case = TRUE)) return("Parallel IGRA (T-SPOT)")
  if (grepl("QFT-GIT|T-SPOT", nm, ignore.case = TRUE)) return("Sequential IGRA")
  return("Symptom screen / CXR")  # TST is a skin test, NOT an IGRA
}

validate_model_correctness <- function(strategies, event_counts_df, icer_table) {
  errors <- character(0)

  # --- 1. IGRA classification invariants ---
  # TST is NOT an IGRA. Only QFT-GIT and T-SPOT.TB are IGRAs.
  all_names <- event_counts_df$strategy

  tst_strategies <- all_names[grepl("TST", all_names, ignore.case = FALSE)]
  for (nm in tst_strategies) {
    grp <- classify_igra(nm)
    if (grp != "Non-IGRA (CXR / symptom screen / TST)")
      errors <- c(errors, sprintf("classify_igra: TST strategy '%s' -> '%s' (should be Non-IGRA)", nm, grp))
    grp2 <- classify_pathway(nm)
    if (grp2 %in% c("Sequential IGRA", "Parallel IGRA (QFT)", "Parallel IGRA (T-SPOT)"))
      errors <- c(errors, sprintf("classify_pathway: TST strategy '%s' -> '%s' (TST is not an IGRA)", nm, grp2))
  }

  # All Parallel strategies must be in the Parallel IGRA group
  parallel_names <- all_names[grepl("^Parallel", all_names)]
  for (nm in parallel_names) {
    grp <- classify_igra(nm)
    if (!grepl("Parallel IGRA", grp))
      errors <- c(errors, sprintf("classify_igra: Parallel strategy '%s' -> '%s' (should be Parallel IGRA)", nm, grp))
    grp2 <- classify_pathway(nm)
    if (!grepl("Parallel IGRA", grp2))
      errors <- c(errors, sprintf("classify_pathway: Parallel strategy '%s' -> '%s' (should be Parallel IGRA)", nm, grp2))
  }

  # Xpert and Ultra (without QFT/T-SPOT/Parallel) must be Non-IGRA
  mol_only <- all_names[grepl("(Xpert|Ultra)$", all_names) &
                         !grepl("QFT|T-SPOT|^Parallel", all_names)]
  for (nm in mol_only) {
    grp <- classify_igra(nm)
    if (grp != "Non-IGRA (CXR / symptom screen / TST)")
      errors <- c(errors, sprintf("classify_igra: Molecular-only strategy '%s' -> '%s' (should be Non-IGRA)", nm, grp))
  }

  # --- 2. Colour palette coverage ---
  # Every named strategy in a palette must exist in the data; no silent grey points
  palettes_to_check <- list(
    ceac_strat_cols    = c("No screening", "Cough+CXR (TB sx)",
                           "Symptom screen+CXR", "Parallel Sx+QFT (Ultra)"),
    ce_ellipse_cols    = c("Cough+CXR (TB sx)", "Symptom screen+CXR",
                           "Parallel Sx+QFT (Ultra)"),
    strat_cols_ov      = c("Parallel Cough+QFT (Ultra)", "Parallel Sx+QFT (Ultra)",
                           "Parallel Cough+T-SPOT (Ultra)", "Parallel Sx+T-SPOT (Ultra)"),
    strat_cols_uptake  = c("Parallel Cough+QFT (Ultra)", "Parallel Sx+QFT (Ultra)",
                           "Parallel Cough+T-SPOT (Ultra)", "Parallel Sx+T-SPOT (Ultra)")
  )
  for (pal_name in names(palettes_to_check)) {
    for (strat in palettes_to_check[[pal_name]]) {
      if (!strat %in% all_names)
        errors <- c(errors, sprintf("Palette '%s' key '%s' not found in strategy names", pal_name, strat))
    }
  }

  # --- 3. Frontier ICERs are cost-effective at NICE threshold ---
  frontier <- icer_table[icer_table$dominance %in% c("ref", "non-dominated"), ]
  for (i in seq_len(nrow(frontier))) {
    icer <- frontier$icer[i]
    if (!is.na(icer) && icer > 25000)
      errors <- c(errors, sprintf("Frontier strategy '%s' has ICER L%s > L25,000/QALY",
                                  frontier$strategy[i], format(round(icer), big.mark=",")))
  }

  # --- 4. No screening must be the lowest-cost strategy ---
  ns_cost <- icer_table$cost[icer_table$strategy == "No screening"]
  if (length(ns_cost) == 0)
    errors <- c(errors, "'No screening' strategy not found in icer_table")
  else if (any(icer_table$cost < ns_cost - 1))
    errors <- c(errors, sprintf("Strategy cheaper than No screening found: %s (L%.0f vs L%.0f)",
                                icer_table$strategy[which.min(icer_table$cost)],
                                min(icer_table$cost), ns_cost))

  # --- 5. Dominance values are valid ---
  valid_dom <- c("ref", "non-dominated", "simply dominated", "extendedly dominated")
  bad_dom <- setdiff(unique(icer_table$dominance), valid_dom)
  if (length(bad_dom) > 0)
    errors <- c(errors, paste("Invalid dominance values:", paste(bad_dom, collapse=", ")))

  # --- 6. Strategy count ---
  n_strat <- nrow(icer_table)
  if (n_strat != 43)
    errors <- c(errors, sprintf("Expected 43 strategies, found %d", n_strat))

  # --- Report ---
  if (length(errors) > 0) {
    cat("\n")
    cat("================================================================================\n")
    cat("  MODEL CORRECTNESS VALIDATION FAILED\n")
    cat("================================================================================\n")
    for (e in errors) cat(sprintf("  x %s\n", e))
    cat("\n")
    stop("Validation failed -- fix errors above before proceeding.", call. = FALSE)
  }
  cat("- Model correctness validation passed (", length(all_names), " strategies checked)\n", sep="")
}

validate_model_correctness(strategies, event_counts_df, icer_table)

ltbi_gap_df <- event_counts_df %>%
  mutate(
    pathway_type = sapply(strategy, classify_pathway),
    pathway_type = factor(pathway_type,
      levels = c("No screening", "Symptom screen / CXR",
                 "Sequential IGRA", "Parallel IGRA (QFT)", "Parallel IGRA (T-SPOT)"))
  ) %>%
  arrange(ltbi_detected) %>%
  mutate(strategy_fct = factor(strategy, levels = strategy))

pathway_pal <- c(
  "No screening"          = "#aaaaaa",
  "Symptom screen / CXR"  = "#e8a0b8",
  "Sequential IGRA"       = "#c090d0",
  "Parallel IGRA (QFT)"   = "#9b5fc0",
  "Parallel IGRA (T-SPOT)"= "#5a1a7a"
)

# Annotation: median sequential IGRA LTBI detected (to label the gap)
seq_igra_median <- ltbi_gap_df %>%
  filter(pathway_type == "Sequential IGRA") %>%
  pull(ltbi_detected) %>% median()

p_ltbi_gap <- ggplot(ltbi_gap_df,
    aes(x = ltbi_detected, y = strategy_fct, fill = pathway_type)) +
  geom_col(width = 0.75, colour = NA) +
  geom_vline(xintercept = seq_igra_median, linetype = "dashed",
             colour = "#c090d0", linewidth = 0.5) +
  annotate("text",
           x = seq_igra_median + 300, y = 3,
           label = sprintf("Sequential\nmax ~%.0f/100k", max(ltbi_gap_df$ltbi_detected[ltbi_gap_df$pathway_type == "Sequential IGRA"])),
           hjust = 0, size = 2.8, colour = "#9b5fc0") +
  annotate("segment",
           x = 1000, xend = 13000, y = 37, yend = 37,
           arrow = arrow(ends = "both", length = unit(0.15, "cm"), type = "open"),
           colour = "#5a1a7a", linewidth = 0.6) +
  annotate("text",
           x = 7000, y = 38.5,
           label = "~147\u00d7 more LTBI detected\n(parallel vs sequential IGRA)",
           hjust = 0.5, size = 2.8, colour = "#5a1a7a", fontface = "bold") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    labels = scales::comma,
    breaks = c(0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000)
  ) +
  scale_fill_manual(values = pathway_pal, name = "Pathway type") +
  theme_minimal(base_size = 10) +
  labs(
    x        = "LTBI cases detected at entry (per 100,000 migrants screened)",
    y        = NULL,
    title    = NULL,
    subtitle = paste0(
      "Sequential strategies (incl. IGRA applied post-symptom-screen): ~72\u2013116/100k\n",
      "Parallel IGRA strategies (universal IGRA): 9,980\u201311,526/100k (specificity-corrected; Pai 2008 BCG)"
    )
  ) +
  theme(
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 9, colour = "grey40"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position    = "bottom",
    legend.text        = element_text(size = 9),
    axis.text.y        = element_text(size = 7.5)
  )

ggsave("output/ltbi_detection_gap.png", p_ltbi_gap,
       width = 11, height = 13, dpi = 300)
cat("Saved: output/ltbi_detection_gap.png\n")

# =============================================================================
# ONE-WAY DETERMINISTIC SENSITIVITY ANALYSIS (DSA)
#
# Each param varied to 2.5th/97.5th percentile (others fixed); NMB at £25k WTP.
# Tornado diagrams for 3 frontier strategies:
# (1) Cough+CXR (TB sx)       seq. ICER £1,179/QALY
# (2) Symptom screen+CXR      seq. ICER £3,822/QALY
# (3) Parallel Sx+QFT (Ultra) seq. ICER £7,745/QALY
#
# Note: only Markov model parameters (costs, QALYs, transition probabilities)
# are varied here. Effective active TB detection rates (eff_active_*) and
# background prevalence values are addressed separately in scenario analyses.
# =============================================================================
# =============================================================================
# IGRA SPECIFICITY SENSITIVITY ANALYSIS
#
# Base case uses test-specific IGRA specificities from Pai et al. 2008 meta-analysis
# (Ann Intern Med 149:177-184; BCG-vaccinated populations):
#   QFT-TB Gold Plus: 0.96 (95% CI 0.94-0.98)   [config: igra_spec_qft]
#   T-SPOT.TB:        0.93 (95% CI 0.86-1.00)    [config: igra_spec_tspt]
#
# True LTBI detection = Sensitivity × LTBI_pool − (1 − Specificity) × Uninfected
# At base case: QFT = 11,526/100k; T-SPOT = 9,980/100k
#
# This SA applies a uniform specificity sweep to both tests to show sensitivity
# across the plausible range:
#   0.90 — conservative lower bound
#   0.93 — T-SPOT.TB base case (Pai 2008 BCG-vaccinated)
#   0.96 — QFT-TB base case    (Pai 2008 BCG-vaccinated)
#   1.00 — naive upper bound (no false positives; shown for comparison only)
#
# For each spec value, 8 parallel IGRA strategies are re-initialised with adjusted
# LTBI detection; non-parallel strategy results held fixed from base case run.
# =============================================================================
cat("\n================================================================================\n")
cat("         IGRA SPECIFICITY SENSITIVITY — PARALLEL LTBI DETECTION                \n")
cat("================================================================================\n\n")

igra_spec_values  <- c(0.90, 0.93, 0.96, 1.00)
n_ltbi_global     <- n_c * prev_ltbi        # 17,800 per 100k
n_uninfected_100k <- n_c * prev_uninfected  # 81,200 per 100k

cat(sprintf("  %-8s  %-20s  %-20s  %-20s\n",
    "Spec.", "QFT LTBI det./100k", "T-SPOT LTBI det./100k", "Note"))
cat(strrep("-", 75), "\n")
for (spec in igra_spec_values) {
  ltbi_qft_s  <- n_ltbi_global * igra_sens_qft  - (1 - spec) * n_uninfected_100k
  ltbi_tspt_s <- n_ltbi_global * igra_sens_tspt - (1 - spec) * n_uninfected_100k
  note <- if (spec == 0.96) "QFT base case (Pai 2008 BCG-vaccinated)" else
          if (spec == 0.93) "T-SPOT base case (Pai 2008 BCG-vaccinated)" else
          if (spec == 1.00) "naive upper bound (no FP; shown for comparison)" else
          "conservative lower bound"
  cat(sprintf("  %-6.2f  %-20.0f  %-20.0f  %s\n",
      spec, ltbi_qft_s, ltbi_tspt_s, note))
}
cat("\n")

# Run frontier ICERs for each specificity level
# (Works directly from strategies list — no dependency on build_strategies internals)
igra_spec_results <- lapply(igra_spec_values, function(spec) {
  ltbi_qft_s  <- max(0, n_ltbi_global * igra_sens_qft  - (1 - spec) * n_uninfected_100k)
  ltbi_tspt_s <- max(0, n_ltbi_global * igra_sens_tspt - (1 - spec) * n_uninfected_100k)

  # Re-create only the 8 parallel IGRA strategies with adjusted LTBI detection
  # Adjust LatentDiagnosed + LatentUndiagnosed in their init vectors directly
  strats_spec <- strategies
  parallel_keys <- names(strats_spec)[sapply(names(strats_spec), function(nm)
    !is.null(strats_spec[[nm]]$prog_cost_type))]
  for (k in parallel_keys) {
    ltbi_s <- if (grepl("tspt", k, ignore.case = TRUE)) ltbi_tspt_s else ltbi_qft_s
    strats_spec[[k]]$init["LatentDiagnosed"]   <- ltbi_s
    strats_spec[[k]]$init["LatentUndiagnosed"] <- n_ltbi_global - ltbi_s
  }

  # Run only affected strategies; copy base results for non-parallel
  spec_results <- strategy_results
  for (k in parallel_keys) {
    r <- run_strategy(strats_spec[[k]])
    idx <- which(sapply(spec_results, function(x) x$strategy_name == r$strategy_name))
    if (length(idx) > 0) spec_results[[idx]] <- r else spec_results <- c(spec_results, list(r))
  }

  icer_s <- calculate_icer(spec_results)
  icer_s$igra_spec <- spec
  icer_s
})

igra_spec_df <- bind_rows(igra_spec_results)
write.csv(igra_spec_df %>% mutate(cost_per_person = cost / n_c, qaly_per_person = qaly / n_c),
          "output/csv/igra_specificity_sensitivity.csv", row.names = FALSE)
cat("Saved: output/csv/igra_specificity_sensitivity.csv\n\n")

# Print frontier comparison across specificity values
frontier_strats_spec <- c("No screening", "Cough+CXR (TB sx)",
                           "Symptom screen+CXR", "Parallel Sx+QFT (Ultra)",
                           "Parallel Sx+T-SPOT (Ultra)")
cat(sprintf("Efficient frontier ICERs by IGRA specificity (QFT base=0.96; T-SPOT base=0.93; Pai 2008 BCG):\n"))
cat(sprintf("  %-32s  %10s  %10s  %10s  %10s\n",
    "Strategy", "Spec=0.90", "Spec=0.93", "Spec=0.96", "Spec=1.00"))
cat(sprintf("  %-32s  %10s  %10s  %10s  %10s\n",
    "", "(conserv.)", "(T-SPOT)", "(QFT base)", "(naive UB)"))
cat(strrep("-", 80), "\n")
for (st in frontier_strats_spec) {
  icers <- sapply(igra_spec_values, function(spec) {
    r <- igra_spec_df[igra_spec_df$strategy == st & abs(igra_spec_df$igra_spec - spec) < 1e-9, ]
    if (nrow(r) == 0 || is.na(r$sequential_icer)) return("—")
    paste0("\u00a3", format(round(r$sequential_icer), big.mark = ","))
  })
  cat(sprintf("  %-32s  %10s  %10s  %10s  %10s\n", st, icers[1], icers[2], icers[3], icers[4]))
}
cat("\n")

# IGRA specificity sensitivity plot
spec_plot_strats <- c("Cough+CXR (TB sx)", "Symptom screen+CXR",
                      "Parallel Sx+QFT (Ultra)", "Parallel Cough+QFT (Ultra)")
spec_pal <- c(
  "Cough+CXR (TB sx)"          = "#c84b6a",
  "Symptom screen+CXR"         = "#9b5fc0",
  "Parallel Sx+QFT (Ultra)"    = "#5a1a7a",
  "Parallel Cough+QFT (Ultra)" = "#8b3a8a"
)
spec_df_plot <- igra_spec_df %>%
  filter(strategy %in% spec_plot_strats, !is.na(sequential_icer)) %>%
  mutate(
    dominated = dominance %in% c("simply dominated", "extendedly dominated"),
    line_type = ifelse(grepl("^Parallel", strategy), "solid", "dashed")
  )

p_igra_spec <- ggplot(spec_df_plot,
                      aes(x = igra_spec, y = sequential_icer, colour = strategy)) +
  geom_line(aes(linetype = line_type), linewidth = 0.9) +
  geom_point(aes(shape = dominated), size = 3.2, alpha = 0.9) +
  scale_colour_manual(values = spec_pal, name = "Strategy") +
  scale_linetype_identity() +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 4),
    labels = c("FALSE" = "On frontier / non-dominated", "TRUE" = "Dominated at this spec."),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = igra_spec_values,
    labels = c("0.90\n(conservative)", "0.93\n(T-SPOT base,\nPai 2008)",
               "0.96\n(QFT base,\nPai 2008)", "1.00\n(naive UB)")
  ) +
  scale_y_continuous(labels = function(x) paste0("\u00a3", format(round(x), big.mark = ","))) +
  theme_minimal(base_size = 13) +
  labs(
    x        = "IGRA specificity (BCG-vaccinated populations; Pai 2008)",
    y        = "Sequential ICER (\u00a3/QALY vs next cheapest non-dominated strategy)",
    subtitle = paste0(
      "Dashed = non-IGRA strategies (unaffected by IGRA specificity) | ",
      "Solid = parallel IGRA strategies\n",
      "X marks = strategy dominated at that specificity value | ",
      "Base case: QFT 0.96, T-SPOT 0.93 (BCG-vaccinated; Pai et al. 2008)"
    )
  ) +
  theme(
    plot.subtitle    = element_text(color = "grey40", size = 9),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    legend.box       = "vertical"
  )

ggsave("output/igra_specificity_sensitivity.png", p_igra_spec,
       width = 11, height = 7, dpi = 300)
cat("Saved: output/igra_specificity_sensitivity.png\n\n")

cat("\n")
cat("================================================================================\n")
cat("             ONE-WAY DETERMINISTIC SENSITIVITY ANALYSIS                         \n")
cat("================================================================================\n\n")

# Find a strategy object by its Excel row name
find_strategy <- function(strats, excel_nm) {
  idx <- which(sapply(strats, function(s) identical(s$excel_name, excel_nm)))
  if (length(idx) == 0) stop(paste("Strategy not found:", excel_nm))
  strats[[idx[1]]]
}

run_owsa <- function(strat_test, strat_ref = NULL) {

  # Restrict DSA to parameters with specified distributions (costs, QALYs,
  # transition probabilities). Fixed structural parameters (n_t, n_c, discount)
  # and detection rates are excluded from one-way variation here.
  dsa_params <- config %>%
    filter(distribution != "fixed", distribution != "", !is.na(distribution)) %>%
    filter(grepl("^(Cost_|qaly_|p_)", parameter))

  if (is.null(strat_ref)) strat_ref <- create_diagnostic_strategies()[["No_screening"]]
  base_result_test <- run_strategy(strat_test)
  base_result_ref  <- run_strategy(strat_ref)
  base_inc_cost    <- base_result_test$total_cost - base_result_ref$total_cost
  base_inc_qaly    <- base_result_test$total_qaly - base_result_ref$total_qaly
  base_icer        <- base_inc_cost / base_inc_qaly
  base_nmb         <- base_inc_qaly * 25000 - base_inc_cost  # NMB at £25,000/QALY WTP (NICE updated Dec 2025)

  # Parallelised one-way SA: each parameter's lo/hi pair is independent,
  # so all rows are computed in parallel with mclapply.
  n_cores_dsa <- max(1L, detectCores() - 1L)

  raw_results <- mclapply(seq_len(nrow(dsa_params)), function(r) {
    pname <- dsa_params$parameter[r]
    pval  <- as.numeric(dsa_params$value[r])
    pdist <- dsa_params$distribution[r]
    p1    <- as.numeric(dsa_params$dist_param1[r])
    p2    <- as.numeric(dsa_params$dist_param2[r])

    if (pdist == "beta") {
      lo <- qbeta(0.025, shape1 = p1, shape2 = p2)
      hi <- qbeta(0.975, shape1 = p1, shape2 = p2)
    } else if (pdist == "gamma") {
      lo <- qgamma(0.025, shape = p1, scale = p2)
      hi <- qgamma(0.975, shape = p1, scale = p2)
    } else {
      return(NULL)  # skip fixed / unsupported distributions
    }

    params_lo          <- paramsData
    params_lo[[pname]] <- lo
    res_test_lo  <- run_strategy(strat_test, params = params_lo)
    res_ref_lo   <- run_strategy(strat_ref,  params = params_lo)
    nmb_lo       <- (res_test_lo$total_qaly - res_ref_lo$total_qaly) * 25000 -
                    (res_test_lo$total_cost  - res_ref_lo$total_cost)

    params_hi          <- paramsData
    params_hi[[pname]] <- hi
    res_test_hi  <- run_strategy(strat_test, params = params_hi)
    res_ref_hi   <- run_strategy(strat_ref,  params = params_hi)
    nmb_hi       <- (res_test_hi$total_qaly - res_ref_hi$total_qaly) * 25000 -
                    (res_test_hi$total_cost  - res_ref_hi$total_cost)

    tibble(
      parameter   = pname,
      param_label = dsa_params$parameter_name[r],
      base_value  = pval,
      lo_value    = lo,
      hi_value    = hi,
      nmb_lo      = nmb_lo,
      nmb_hi      = nmb_hi,
      nmb_base    = base_nmb,
      nmb_range   = abs(nmb_hi - nmb_lo)
    )
  }, mc.cores = n_cores_dsa)

  dsa_df <- bind_rows(Filter(Negate(is.null), raw_results)) %>% arrange(desc(nmb_range))
  return(list(dsa_df = dsa_df, base_nmb = base_nmb, base_icer = base_icer))
}

# Helper: build and save a tornado diagram from a run_owsa() result.
# bar_col  = fill colour for bars; cap_col = endpoint marker colour.
# Short human-readable labels for tornado diagram y-axis
dsa_short_labels <- c(
  "QALY weight - Latent TB undiagnosed"              = "QALY: latent undiagnosed",
  "QALY weight - Uninfected"                         = "QALY: uninfected",
  "QALY weight - Latent TB not treated"              = "QALY: latent untreated",
  "Prob: Latent undiagnosed to Active undiagnosed"   = "Reactivation rate",
  "Prob: Active undiagnosed to Dead"                 = "TB case fatality (undiag.)",
  "Prob: Active undiagnosed to diagnosed"            = "Spontaneous TB diagnosis",
  "QALY weight - Latent TB lost to follow-up"        = "QALY: latent LTFU",
  "Prob: Latent diagnosed to not treated"            = "LTBI treatment non-initiation",
  "Prob: Latent not treated to Active undiagnosed"   = "Reactivation (untreated LTBI)",
  "Prob: Active treatment to completed"              = "Active TB treatment completion",
  "Prob: Active treatment to discontinued"           = "Active TB treatment discontinuation",
  "QALY weight - Latent TB treatment completed"      = "QALY: LTBI treated (completed)",
  "Prob: Latent diagnosed to treatment"              = "LTBI treatment initiation",
  "QALY weight - Latent TB diagnosed"                = "QALY: LTBI diagnosed",
  "Prob: Latent undiagnosed to diagnosed"            = "Incidental LTBI diagnosis",
  "QALY weight - Latent TB treatment discontinued"   = "QALY: LTBI tx discontinued",
  "Prob: Latent treatment to lost to follow-up"      = "LTBI treatment LTFU",
  "QALY weight - Active TB undiagnosed"              = "QALY: active TB (undiag.)",
  "Prob: Latent LTFU to Active undiagnosed"          = "Reactivation (LTBI LTFU)",
  "Prob: Active discontinued to undiagnosed"         = "Active TB: discont. → undiag.",
  "Prob: Active diagnosed to not treated"            = "Active TB: diagnosed → untreated",
  "Prob: Active not treated to diagnosed"            = "Active TB: untreated → diagnosed",
  "Prob: Active treated to Dead"                     = "Active TB mortality (on tx)",
  "Prob: Active treatment to lost to follow-up"      = "Active TB tx LTFU",
  "Prob: Latent treatment to discontinued"           = "LTBI tx discontinuation",
  "Prob: Active discontinued to Dead"                = "Active TB mortality (discont.)",
  "Cost - Active TB under treatment"                 = "Cost: active TB treatment",
  "Prob: Active not treated to Dead"                 = "Active TB mortality (untreated)",
  "Prob: Active diagnosed to Dead"                   = "Active TB mortality (diagnosed)",
  "Prob: Active diagnosed to treatment"              = "Active TB treatment initiation",
  "QALY weight - Active TB under treatment"          = "QALY: active TB (on tx)",
  "QALY weight - Active TB treatment completed"      = "QALY: active TB (completed)",
  "QALY weight - Active TB treatment discontinued"   = "QALY: active TB (discont.)",
  "QALY weight - Active TB diagnosed"                = "QALY: active TB (diagnosed)",
  "QALY weight - Active TB not treated"              = "QALY: active TB (untreated)",
  "QALY weight - Active TB lost to follow-up"        = "QALY: active TB (LTFU)"
)

plot_tornado_diagram <- function(owsa_result, subtitle_str, bar_col, cap_col, out_file,
                                 top_n_override = NULL, base_size = 13, save = TRUE,
                                 label_width = 40, show_caption = TRUE, show_title = TRUE) {
  top_n  <- if (!is.null(top_n_override)) min(top_n_override, nrow(owsa_result$dsa_df))
             else min(15, nrow(owsa_result$dsa_df))
  t_df   <- owsa_result$dsa_df %>%
    slice_head(n = top_n) %>%
    mutate(
      param_short = dplyr::coalesce(dsa_short_labels[param_label], param_label),
      param_short = ifelse(nchar(param_short) > label_width,
                           paste0(substr(param_short, 1, label_width - 3), "..."),
                           param_short),
      param_short = make.unique(param_short, sep = " "),
      param_short = factor(param_short, levels = rev(param_short)),
      nmb_lo_dev  = nmb_lo - owsa_result$base_nmb,
      nmb_hi_dev  = nmb_hi - owsa_result$base_nmb,
      bar_left    = pmin(nmb_lo_dev, nmb_hi_dev),
      bar_right   = pmax(nmb_lo_dev, nmb_hi_dev)
    )
  p <- ggplot(t_df) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
    geom_segment(aes(y = param_short, yend = param_short,
                     x = bar_left, xend = bar_right),
                 linewidth = 8, color = bar_col, alpha = 0.8) +
    geom_point(aes(y = param_short, x = bar_left),  shape = "|", size = 3, color = cap_col) +
    geom_point(aes(y = param_short, x = bar_right), shape = "|", size = 3, color = cap_col) +
    scale_x_continuous(
      labels = function(x) {
        sign_str <- ifelse(x < 0, "\u2212\u00a3", "\u00a3")
        ax <- abs(x)
        ifelse(ax >= 1e6, paste0(sign_str, round(ax / 1e6, 1), "M"),
        ifelse(ax >= 1e3, paste0(sign_str, round(ax / 1e3, 0), "k"),
                          paste0(sign_str, round(ax))))
      },
      n.breaks = 7
    ) +
    theme_minimal(base_size = base_size) +
    labs(
      x        = "Change in NMB from base case (2023/24 \u00a3)",
      y        = "",
      title    = if (show_title) "One-way deterministic sensitivity analysis" else NULL,
      subtitle = subtitle_str,
      caption  = if (show_caption) "NMB = net monetary benefit; WTP = willingness-to-pay threshold (NICE \u00a325,000/QALY); DSA = deterministic sensitivity analysis.\nAll costs in 2023/24 GBP. Parameters varied individually across 95% plausible range; all others held at base-case values." else NULL
    ) +
    theme(
      plot.title         = element_text(face = "bold", size = base_size + 2),
      plot.subtitle      = element_text(color = "grey40", size = base_size - 1),
      plot.caption       = element_text(color = "grey40", size = max(base_size - 3, 7),
                                        hjust = 0),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y        = element_text(size = base_size - 1)
    )
  if (save) {
    ggsave(out_file, p, width = 14, height = 8, dpi = 300)
    cat(sprintf("Tornado diagram saved to %s\n", out_file))
  }
  invisible(p)
}

cat("Running one-way DSA (Parallel Sx+QFT (Ultra) vs No Screening)...\n")
owsa <- run_owsa(
  strat_test = find_strategy(strategies, "parallel allsx qft_ultra"),
  strat_ref  = strategies[["No_screening"]]
)
dsa_df <- owsa$dsa_df

# Report the 15 parameters with the largest impact on NMB
cat("\nTop 15 most influential parameters on NMB (£25k WTP):\n")
cat(sprintf("Base case NMB: £%s | Base case ICER: £%s\n\n",
            format(round(owsa$base_nmb), big.mark = ","),
            format(round(owsa$base_icer), big.mark = ",")))
cat(sprintf("%-45s %10s %12s %12s %12s\n",
            "Parameter", "Base", "NMB (Low)", "NMB (High)", "Range"))
cat(sprintf("%-45s %10s %12s %12s %12s\n",
            "---------------------------------------------", "----------",
            "------------", "------------", "------------"))
for (i in 1:min(15, nrow(dsa_df))) {
  cat(sprintf("%-45s %10.4f %11s %11s %11s\n",
              substr(dsa_df$param_label[i], 1, 45),
              dsa_df$base_value[i],
              paste0("£", format(round(dsa_df$nmb_lo[i]), big.mark = ",")),
              paste0("£", format(round(dsa_df$nmb_hi[i]), big.mark = ",")),
              paste0("£", format(round(dsa_df$nmb_range[i]), big.mark = ","))))
}

# Save full DSA results table for supplementary materials
write.csv(dsa_df, "output/csv/dsa_results.csv", row.names = FALSE)

# -------------------- Tornado diagram -----------------------------------------
# Displays the top 15 parameters ranked by their impact on NMB. Each bar shows
# how NMB changes when the parameter moves from its lower to upper bound.
# The width of each bar reflects the parameter's contribution to overall
# model uncertainty. Bars are expressed as deviation from the base case NMB
# so that the central reference line (x = 0) represents no change.
p_tornado_sx_qft <- plot_tornado_diagram(owsa,
  subtitle_str = "Parallel Sx+QFT (Ultra) vs No Screening | NMB at \u00a325,000/QALY WTP | Top 15 parameters",
  bar_col = "#9b5fc0", cap_col = "#6a3090",
  out_file = "output/tornado_diagram.png")

# -------------------- Second tornado: Symptom screen+CXR vs No Screening
# Symptom screen+CXR is the second strategy on the efficient frontier
# (seq. ICER £3,822/QALY). This tornado characterises its structural uncertainty.
cat("Running one-way DSA (Symptom screen+CXR vs No Screening)...\n")
owsa_cxr <- run_owsa(
  strat_test = find_strategy(strategies, "anysx_cxr(any)"),
  strat_ref  = strategies[["No_screening"]]
)

cat(sprintf("Base case NMB (Symptom screen+CXR): £%s | Base case ICER: £%s\n\n",
            format(round(owsa_cxr$base_nmb), big.mark = ","),
            format(round(owsa_cxr$base_icer), big.mark = ",")))

p_tornado_cough_qft <- plot_tornado_diagram(owsa_cxr,
  subtitle_str = "Symptom screen+CXR vs No Screening | NMB at \u00a325,000/QALY WTP | Top 15 parameters",
  bar_col = "#c84b6a", cap_col = "#a03060",
  out_file = "output/tornado_diagram_symscrCXR.png")

# -------------------- Third tornado: Cough+CXR (TB sx) vs No Screening --------
# Cough+CXR (TB sx) is the lowest-cost strategy on the efficient frontier
# (seq. ICER ~£1,210/QALY vs No Screening; well within NICE £25,000 threshold).
# This tornado characterises the structural uncertainty around the cheapest
# CE strategy — the non-IGRA entry-level recommendation.
cat("Running one-way DSA (Cough+CXR (TB sx) vs No Screening)...\n")
owsa_cough_cxr <- run_owsa(
  strat_test = find_strategy(strategies, "anycough_cxr(TB)"),
  strat_ref  = strategies[["No_screening"]]
)

cat(sprintf("Base case NMB (Cough+CXR TB sx): £%s | Base case ICER: £%s\n\n",
            format(round(owsa_cough_cxr$base_nmb), big.mark = ","),
            format(round(owsa_cough_cxr$base_icer), big.mark = ",")))

p_tornado_cough_cxr <- plot_tornado_diagram(owsa_cough_cxr,
  subtitle_str = "Cough+CXR (TB sx) vs No Screening | NMB at \u00a325,000/QALY WTP | Top 15 parameters",
  bar_col = "#e06888", cap_col = "#a03060",
  out_file = "output/tornado_diagram_coughCXR_TBsx.png")

# Save DSA results for Cough+CXR and Symptom screen+CXR (for panel regeneration)
write.csv(owsa_cough_cxr$dsa_df, "output/csv/dsa_results_coughcxr.csv",  row.names = FALSE)
write.csv(owsa_cxr$dsa_df,       "output/csv/dsa_results_symscrCXR.csv",  row.names = FALSE)

# -------------------- Tornado panel: all 3 frontier strategies ----------------
# Three panels, one per non-reference frontier strategy, ordered by ascending cost:
#   A. Cough+CXR (TB sx)        — ICER £1,179/QALY
#   B. Symptom screen+CXR       — ICER £3,822/QALY
#   C. Parallel Sx+QFT (Ultra)  — ICER £7,745/QALY
# Top 7 parameters each (top 8 too cramped at 3-panel width); shared caption.
p_tA <- plot_tornado_diagram(owsa_cough_cxr,
  subtitle_str = "A. Cough+CXR (TB sx)",
  bar_col = "#e06888", cap_col = "#a03060",
  top_n_override = 7, base_size = 9, label_width = 28, save = FALSE, show_caption = FALSE, show_title = FALSE)
p_tB <- plot_tornado_diagram(owsa_cxr,
  subtitle_str = "B. Symptom screen+CXR",
  bar_col = "#c84b6a", cap_col = "#a03060",
  top_n_override = 7, base_size = 9, label_width = 28, save = FALSE, show_caption = FALSE, show_title = FALSE)
p_tC <- plot_tornado_diagram(owsa,
  subtitle_str = "C. Parallel Sx+QFT (Ultra)",
  bar_col = "#9b5fc0", cap_col = "#6a3090",
  top_n_override = 7, base_size = 9, label_width = 28, save = FALSE, show_caption = FALSE, show_title = FALSE)

p_tornado_panel <- (p_tA | p_tB | p_tC) +
  patchwork::plot_annotation(
    title   = NULL,
    caption = "NMB at \u00a325,000/QALY WTP threshold vs No Screening. Bars show change in NMB when each parameter varies between its lower and upper bound; all other parameters held at base-case values. Top 7 parameters by NMB impact shown per strategy.",
    theme   = theme(plot.title   = element_text(face = "bold", size = 14),
                    plot.caption = element_text(color = "grey40", size = 9))
  )
ggsave("output/tornado_panel.png", p_tornado_panel, width = 15, height = 6, dpi = 300)
cat("Tornado panel saved to output/tornado_panel.png\n")

# =============================================================================
# IGRA TESTING PROGRAMME COST SENSITIVITY ANALYSIS — PARALLEL IGRA STRATEGIES
#
# Analytical one-way sensitivity on IGRA testing programme delivery cost.
# Programme delivery cost is a one-time lump-sum test cost that does not affect
# Markov transitions, so ICER changes linearly with programme cost:
#
#   ICER(pc) = [base_inc_cost + (pc − pc_base)] / base_inc_qaly
#
# where base_inc_cost and base_inc_qaly are cohort-level totals vs No Screening.
# No model re-runs are required.
#
# Range swept: £1–£10 per person (= £100,000–£1,000,000 per 100,000 migrants),
# covering the reviewer-requested £150k–£800k range with margin on each side.
# Base case points: cough strategies £3/pp (£300k), allsx strategies £5/pp (£500k).
# =============================================================================
cat("\n================================================================================\n")
cat("    IGRA TESTING PROGRAMME COST SENSITIVITY — PARALLEL IGRA STRATEGIES          \n")
cat("================================================================================\n\n")

parallel_keys_all <- names(strategies)[sapply(names(strategies), function(nm)
  !is.null(strategies[[nm]]$prog_cost_type))]

ref_sr       <- strategy_results[[which(sapply(strategy_results, function(x) x$strategy_name == "No screening"))]]
ref_cost_tot <- ref_sr$total_cost
ref_qaly_tot <- ref_sr$total_qaly

prog_cost_pp_seq <- seq(0.5, 60, by = 0.5)  # £0.50–£60 per person (£50k–£6M per 100k)
# Extends beyond reviewer-requested £1.50–£8.00/pp to show that strategies remain
# cost-effective even at 7× the base case programme cost (breakeven ≈ £1,300/pp; far off-chart)

prog_cost_icer_df <- bind_rows(lapply(parallel_keys_all, function(k) {
  strat  <- strategies[[k]]
  sr_idx <- which(sapply(strategy_results, function(x) x$strategy_name == strat$name))
  if (length(sr_idx) == 0) return(NULL)
  sr <- strategy_results[[sr_idx]]

  pc_base      <- strat$prog_cost_total         # base programme delivery cost in £ total
  base_ic_cost <- sr$total_cost - ref_cost_tot  # incremental cost vs No Screening at base programme cost
  base_ic_qaly <- sr$total_qaly - ref_qaly_tot  # incremental QALY (fixed — programme cost doesn't affect QALY)
  if (base_ic_qaly <= 0) return(NULL)           # skip dominated at base

  bind_rows(lapply(prog_cost_pp_seq, function(pc_pp) {
    pc          <- pc_pp * n_c               # total programme delivery cost £ for this scenario
    new_ic_cost <- base_ic_cost + (pc - pc_base)
    tibble(
      strategy        = strat$name,
      prog_cost_type  = strat$prog_cost_type,
      prog_cost_pp    = pc_pp,
      prog_cost_total = pc,
      icer            = new_ic_cost / base_ic_qaly
    )
  }))
}))

write.csv(prog_cost_icer_df, "output/csv/igra_programme_cost_sensitivity.csv", row.names = FALSE)
cat("Saved: output/csv/igra_programme_cost_sensitivity.csv\n")

# Print breakeven programme cost for each strategy (ICER = £25,000/QALY)
cat("\nBreakeven programme delivery cost (ICER = £25,000/QALY) for parallel IGRA strategies:\n")
cat(sprintf("  %-34s  %12s  %12s  %12s\n", "Strategy", "ICER (base)", "Breakeven £/pp", "Breakeven £/100k"))
cat(strrep("-", 80), "\n")
for (k in parallel_keys_all) {
  strat  <- strategies[[k]]
  df_k   <- prog_cost_icer_df %>% filter(strategy == strat$name)
  if (nrow(df_k) == 0) next
  sr_idx <- which(sapply(strategy_results, function(x) x$strategy_name == strat$name))
  sr     <- strategy_results[[sr_idx]]
  base_icer_k    <- (sr$total_cost - ref_cost_tot) / (sr$total_qaly - ref_qaly_tot)
  base_ic_qaly_k <- sr$total_qaly - ref_qaly_tot
  # Breakeven: ICER(pc) = 25000 => pc = pc_base + (25000 - base_icer_k) * base_ic_qaly
  breakeven_pc   <- strat$prog_cost_total + (25000 - base_icer_k) * base_ic_qaly_k
  cat(sprintf("  %-34s  %11s  %12.2f  %12.0f\n",
              strat$name,
              paste0("\u00a3", format(round(base_icer_k), big.mark = ",")),
              breakeven_pc / n_c,
              breakeven_pc))
}
cat("\n")

# --- Plot: ICER vs programme delivery cost per person ---
# Focus on Ultra confirmatory strategies (the frontier ones)
plot_strats <- c("Parallel Cough+QFT (Ultra)", "Parallel Sx+QFT (Ultra)",
                 "Parallel Cough+T-SPOT (Ultra)", "Parallel Sx+T-SPOT (Ultra)")

strat_cols_ov <- c(
  "Parallel Cough+QFT (Ultra)"    = "#9b5fc0",
  "Parallel Sx+QFT (Ultra)"       = "#5a1a7a",
  "Parallel Cough+T-SPOT (Ultra)" = "#d06090",
  "Parallel Sx+T-SPOT (Ultra)"    = "#c84b6a"
)

prog_cost_plot_df <- prog_cost_icer_df %>%
  filter(strategy %in% plot_strats) %>%
  mutate(strategy = factor(strategy, levels = plot_strats))

# Base case annotation points (cough = £3/pp, allsx = £5/pp)
base_pts <- prog_cost_plot_df %>%
  mutate(is_base = (prog_cost_type == "cough" & abs(prog_cost_pp - 3) < 0.1) |
                   (prog_cost_type == "allsx" & abs(prog_cost_pp - 5) < 0.1)) %>%
  filter(is_base)

# Reviewer-requested range shading: £1.50–£8.00/pp
reviewer_lo <- 1.5
reviewer_hi <- 8.0

p_prog_cost <- ggplot(prog_cost_plot_df, aes(x = prog_cost_pp, y = icer,
                                              colour = strategy, group = strategy)) +
  # Shade reviewer-requested range
  annotate("rect", xmin = reviewer_lo, xmax = reviewer_hi, ymin = -Inf, ymax = Inf,
           fill = "#f5e8f5", alpha = 0.5) +
  annotate("text", x = (reviewer_lo + reviewer_hi) / 2, y = 28500,
           label = "Reviewer-requested\nrange (£1.50\u2013£8.00/pp)", size = 2.6,
           colour = "#9070a0", hjust = 0.5) +
  geom_line(linewidth = 1.1) +
  # Base case vertical lines
  geom_vline(xintercept = 3, linetype = "dotdash", colour = "#7a50a8", linewidth = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 5, linetype = "dotdash", colour = "#c84b6a", linewidth = 0.5, alpha = 0.8) +
  # WTP thresholds
  geom_hline(yintercept = 25000, linetype = "dashed", colour = "grey30", linewidth = 0.8) +
  geom_hline(yintercept = 35000, linetype = "dotted", colour = "grey50", linewidth = 0.7) +
  annotate("text", x = 59, y = 26200, label = "\u00a325,000/QALY (NICE lower)", hjust = 1, size = 3.0, colour = "grey30") +
  annotate("text", x = 59, y = 36200, label = "\u00a335,000/QALY (NICE upper)", hjust = 1, size = 3.0, colour = "grey50") +
  # Base case labels
  annotate("text", x = 3.5, y = 1200, label = "Base (cough)\n£3/pp", size = 2.6,
           colour = "#7a50a8", hjust = 0) +
  annotate("text", x = 5.5, y = 1200, label = "Base (all-sx)\n£5/pp", size = 2.6,
           colour = "#c84b6a", hjust = 0) +
  # Breakeven off-chart annotation
  annotate("text", x = 55, y = 3500,
           label = "Breakeven \u2248 £1,300/pp\n(off-chart; all strategies\nremain CE far beyond\nplausible programme cost range)",
           size = 2.5, colour = "grey40", hjust = 1) +
  # Base case dots
  geom_point(data = base_pts, aes(colour = strategy), size = 3, shape = 21,
             fill = "white", stroke = 1.5) +
  scale_colour_manual(values = strat_cols_ov, name = "Strategy") +
  scale_x_continuous(
    breaks = c(0.5, seq(5, 60, by = 5)),
    labels = function(x) paste0("\u00a3", x, "\n(", format(round(x * 100), big.mark = ","), "k/100k)")
  ) +
  scale_y_continuous(
    labels  = function(y) paste0("\u00a3", format(round(y), big.mark = ",")),
    limits  = c(0, 40000),
    breaks  = seq(0, 40000, by = 5000)
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title    = NULL,
    subtitle = paste0("Parallel IGRA (Ultra) strategies vs No Screening | Analytical one-way SA\n",
                      "Purple shading = reviewer-requested range | Dashed = NICE WTP thresholds | ",
                      "Dots = base case programme cost (cough: \u00a33/pp; all-sx: \u00a35/pp)"),
    x = "IGRA testing programme delivery cost per person [\u00a3/pp] (equivalent per 100,000 migrants)",
    y = "Incremental ICER vs No Screening (\u00a3/QALY)"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.minor = element_blank()
  ) +
  guides(colour = guide_legend(nrow = 2))

ggsave("output/igra_programme_cost_sensitivity.png", p_prog_cost, width = 12, height = 7, dpi = 300)
cat("Saved: output/igra_programme_cost_sensitivity.png\n")

cat("\n================================================================================\n")

# =============================================================================
# IGRA UPTAKE SENSITIVITY (PARALLEL STRATEGIES)
#
# Base case assumes 100% of LTBI-positive migrants accept the IGRA test
# (LatentDiagnosed = n_ltbi × IGRA sensitivity). In practice uptake in UK
# migrant screening programmes is typically 60-80% (UKHSA 2024).
#
# This sensitivity sweeps uptake from 50% to 100% in 5% steps.
# For each level, LatentDiagnosed is scaled proportionally and
# LatentUndiagnosed compensated. Costs are held fixed (programme delivery cost is a
# programme-level cost independent of individual uptake). Pairwise ICERs
# vs No Screening are recomputed for all 8 parallel IGRA strategies.
# =============================================================================
cat("\n================================================================================\n")
cat("         IGRA UPTAKE SENSITIVITY (PARALLEL STRATEGIES)                          \n")
cat("================================================================================\n\n")

# Identify parallel IGRA strategies by presence of prog_cost_type field
parallel_keys <- names(strategies)[sapply(names(strategies), function(nm)
  !is.null(strategies[[nm]]$prog_cost_type))]

igra_uptake_vals <- seq(0.50, 1.00, by = 0.05)

# Run No Screening once as reference
res_no_screen_uptake <- run_strategy(strategies[["No_screening"]])

cat(sprintf("Sweeping IGRA uptake: %.0f%% to %.0f%% (%d steps) across %d parallel strategies\n\n",
            min(igra_uptake_vals)*100, max(igra_uptake_vals)*100,
            length(igra_uptake_vals), length(parallel_keys)))

uptake_rows <- lapply(igra_uptake_vals, function(u) {
  strat_results <- mclapply(parallel_keys, function(nm) {
    s <- strategies[[nm]]
    orig_ld <- s$init["LatentDiagnosed"]
    orig_lu <- s$init["LatentUndiagnosed"]
    # Scale diagnosed down by uptake; redirect non-uptakers to undiagnosed
    # (conserves total LTBI count without needing n_ltbi in scope)
    s$init["LatentDiagnosed"]   <- orig_ld * u
    s$init["LatentUndiagnosed"] <- orig_lu + orig_ld * (1 - u)
    # NOTE: test_cost is held constant across uptake levels (conservative assumption).
    # test_cost = sequential_base + programme_cost, where programme_cost (£300k–£500k/100k)
    # is the dominant component and represents fixed programme infrastructure (staffing,
    # logistics, equipment). At lower uptake the IGRA reagent variable cost is overstated,
    # so ICERs are slightly inflated — this provides a pessimistic bound on cost-effectiveness.
    run_strategy(s)
  }, mc.cores = max(1L, detectCores() - 1L))

  lapply(strat_results, function(res) {
    inc_cost <- res$total_cost  - res_no_screen_uptake$total_cost
    inc_qaly <- res$total_qaly  - res_no_screen_uptake$total_qaly
    tibble(
      strategy = res$strategy_name,
      uptake   = u,
      icer_vs_noscreening = if (inc_qaly > 0) inc_cost / inc_qaly else NA_real_,
      cost_pp  = res$total_cost  / n_c,
      qaly_pp  = res$total_qaly  / n_c
    )
  }) %>% bind_rows()
}) %>% bind_rows()

# Print summary at base-case uptake levels (70%, 85%, 100%)
for (u_show in c(0.70, 0.85, 1.00)) {
  cat(sprintf("--- Uptake = %.0f%% ---\n", u_show * 100))
  sub <- uptake_rows %>%
    filter(abs(uptake - u_show) < 0.001) %>%
    arrange(icer_vs_noscreening)
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("  %-38s  ICER vs No Screening: %s\n",
                sub$strategy[i],
                if (is.na(sub$icer_vs_noscreening[i])) "N/A"
                else paste0("\u00a3", format(round(sub$icer_vs_noscreening[i]), big.mark = ","))))
  }
  cat("\n")
}

write.csv(uptake_rows, "output/csv/igra_uptake_sensitivity.csv", row.names = FALSE)
cat("Saved: output/csv/igra_uptake_sensitivity.csv\n\n")

# Plot: ICER vs uptake for the 4 Ultra parallel strategies (most policy-relevant)
ultra_strategies <- uptake_rows %>%
  filter(grepl("Ultra", strategy)) %>%
  mutate(uptake_pct = uptake * 100)

strat_cols_uptake <- c(
  "Parallel Cough+QFT (Ultra)"    = "#9b5fc0",
  "Parallel Sx+QFT (Ultra)"       = "#5a1a7a",
  "Parallel Cough+T-SPOT (Ultra)" = "#d06090",
  "Parallel Sx+T-SPOT (Ultra)"    = "#a03060"
)

p_uptake <- ggplot(ultra_strategies,
                   aes(x = uptake_pct, y = icer_vs_noscreening,
                       colour = strategy, group = strategy)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 100, linetype = "dotdash", colour = "grey50",
             linewidth = 0.5, alpha = 0.8) +
  annotate("text", x = 100, y = 3000,
           label = "Base case\n(100% uptake)", hjust = 1.05, size = 2.8,
           colour = "grey40") +
  geom_hline(yintercept = 25000, linetype = "dashed", colour = "grey30",
             linewidth = 0.8) +
  geom_hline(yintercept = 35000, linetype = "dotted", colour = "grey50",
             linewidth = 0.7) +
  annotate("text", x = 51, y = 26200,
           label = "\u00a325,000/QALY (NICE lower)", hjust = 0, size = 3.0,
           colour = "grey30") +
  annotate("text", x = 51, y = 36200,
           label = "\u00a335,000/QALY (NICE upper)", hjust = 0, size = 3.0,
           colour = "grey50") +
  scale_colour_manual(values = strat_cols_uptake, name = "Strategy") +
  scale_x_continuous(breaks = seq(50, 100, by = 10),
                     labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels  = function(y) paste0("\u00a3", format(round(y), big.mark = ",")),
                     limits  = c(0, 40000),
                     breaks  = seq(0, 40000, by = 5000)) +
  theme_minimal(base_size = 12) +
  labs(
    title    = NULL,
    subtitle = paste0("Pairwise ICER vs No Screening | Parallel IGRA (Ultra) strategies\n",
                      "Base case assumes 100% uptake; UKHSA 2024 real-world range ~60\u201380%"),
    x = "IGRA test uptake (%)",
    y = "ICER vs No Screening (\u00a3/QALY)"
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(colour = guide_legend(nrow = 2))

ggsave("output/igra_uptake_sensitivity.png", p_uptake,
       width = 10, height = 7, dpi = 300)
cat("Saved: output/igra_uptake_sensitivity.png\n")

cat("\n================================================================================\n")

# =============================================================================
# LTBI TREATMENT COMPLETION SENSITIVITY
#
# Base case: p_LatentTreated_LatentCompleted = 0.253/mo
#   Source: Surey et al. 2021 (3HR arm, n=25, London UK): 76% completion / 3 months
#   With p_LTFU = 0.057/mo (also Surey 2021: 4/25 LFU over 12wk) and
#   p_discontinued = 0.02/mo, conditional completion = 0.253/(0.253+0.02+0.057) = 76.7%
#   → consistent with 76% trial primary outcome. Base case is now self-calibrating.
#
# Real-world lower bound: 0.185/mo
#   Source: UKHSA TB in England 2025 (Prevention: England 2024): 55.5% completion
#   in contacts / 3 months = 0.185/mo. Tests sensitivity to lower real-world adherence.
#   Note: UKHSA also reports 11.7% overall completion but this figure has large
#   missing data and is not suitable as a SA bound (see config.csv row 48 note).
#
# Direction of effect: lower completion → fewer people reach LatentCompleted →
#   slightly less QALY benefit from LTBI treatment → frontier ICERs increase.
#   Expected to remain <<£25k.
# =============================================================================
cat("\n================================================================================\n")
cat("    LTBI TREATMENT COMPLETION \u2014 TRIAL vs REAL-WORLD SENSITIVITY               \n")
cat("================================================================================\n\n")

p_ltbi_completion_rw <- as.numeric(config_list["p_ltbi_completion_rw"]) / 3

p_ltfu_base <- as.numeric(config_list["p_LatentTreated_LatentLtfu"])
p_disc_base <- as.numeric(config_list["p_LatentTreated_LatentDiscontinued"])
p_comp_base <- as.numeric(config_list["p_LatentTreated_LatentCompleted"])
cond_completion_base <- 100 * p_comp_base / (p_comp_base + p_disc_base + p_ltfu_base)

cat(sprintf("Base case:  p_LatentTreated_LatentCompleted = %.4f/mo  (Surey 2021: 76%% / 3mo)\n",
            p_comp_base))
cat(sprintf("            p_LatentTreated_LatentLtfu       = %.4f/mo  (Surey 2021: 4/25 LFU 3HR arm)\n",
            p_ltfu_base))
cat(sprintf("            Conditional completion in model  = %.1f%%  (target: 76%% trial rate)\n\n",
            cond_completion_base))
cat(sprintf("Sensitivity: p_LatentTreated_LatentCompleted = %.4f/mo  (UKHSA 2025: 55.5%% / 3mo)\n\n",
            p_ltbi_completion_rw))

params_ltbi_rw <- paramsData
params_ltbi_rw[["p_LatentTreated_LatentCompleted"]] <- p_ltbi_completion_rw

strategy_results_ltbi_rw <- mclapply(strategies, run_strategy,
                                      params   = params_ltbi_rw,
                                      mc.cores = max(1L, detectCores() - 1L))
icer_table_ltbi_rw <- calculate_icer(strategy_results_ltbi_rw)

cat("ICER TABLE \u2014 Sensitivity: real-world LTBI completion (55.5%; base case: 76%):\n")
cat(sprintf("%-22s %11s %10s %11s %14s %15s  %s\n",
            "Strategy", "Cost/person", "QALYs/pp", "Inc.Cost",
            "ICER(vs ref)", "Sequential ICER", "Dominance"))
cat(strrep("-", 108), "\n")
for (i in 1:nrow(icer_table_ltbi_rw)) {
  icer_str     <- if (is.na(icer_table_ltbi_rw$icer[i])) "Ref"
                  else paste0("\u00a3", format(round(icer_table_ltbi_rw$icer[i]), big.mark = ","))
  seq_icer_str <- if (is.na(icer_table_ltbi_rw$sequential_icer[i])) "-"
                  else paste0("\u00a3", format(round(icer_table_ltbi_rw$sequential_icer[i]), big.mark = ","))
  cat(sprintf("%-22s %10s %10.4f %10s %14s %15s  %s\n",
              icer_table_ltbi_rw$strategy[i],
              paste0("\u00a3", format(round(icer_table_ltbi_rw$cost_per_person[i]), big.mark = ",")),
              icer_table_ltbi_rw$qaly_per_person[i],
              paste0("\u00a3", format(round(icer_table_ltbi_rw$inc_cost[i] / n_c), big.mark = ",")),
              icer_str,
              seq_icer_str,
              icer_table_ltbi_rw$dominance[i]))
}
cat("\n")

# Focused frontier comparison (base vs real-world vs calibrated) — 3 scenarios
frontier_names_ltbi <- icer_table %>%
  filter(dominance %in% c("ref", "non-dominated")) %>%
  arrange(cost) %>%
  pull(strategy)

cat(sprintf("FRONTIER COMPARISON \u2014 LTBI Completion: Base (~%.0f%% conditional) vs Real-World (55.5%% UKHSA 2025):\n",
            cond_completion_base))
cat(sprintf("  %-32s  %26s  %26s  %12s\n",
    "Strategy",
    sprintf("Base (~%.0f%%/0.253/mo)", cond_completion_base),
    "RW (55.5%/0.185/mo)",
    "Delta"))
cat(strrep("-", 104), "\n")
for (st in frontier_names_ltbi) {
  r_base <- icer_table[icer_table$strategy == st, ]
  r_rw   <- icer_table_ltbi_rw[icer_table_ltbi_rw$strategy == st, ]
  icer_b <- if (is.na(r_base$sequential_icer)) "Ref" else
            paste0("\u00a3", format(round(r_base$sequential_icer), big.mark = ","), "/QALY")
  icer_r <- if (is.na(r_rw$sequential_icer)) "Ref" else
            paste0("\u00a3", format(round(r_rw$sequential_icer), big.mark = ","), "/QALY")
  delta  <- if (is.na(r_base$sequential_icer) || is.na(r_rw$sequential_icer)) "\u2014" else
            paste0(ifelse(r_rw$sequential_icer - r_base$sequential_icer >= 0, "+", ""),
                   "\u00a3", format(round(r_rw$sequential_icer - r_base$sequential_icer),
                                    big.mark = ","))
  cat(sprintf("  %-32s  %26s  %26s  %12s\n", st, icer_b, icer_r, delta))
}
cat("\n================================================================================\n\n")

# Save combined CSV (2 scenarios: base trial rate, real-world UKHSA 2025)
icer_table_ltbi_combined <- bind_rows(
  icer_table %>%
    mutate(scenario = "trial_76pct",
           cost_per_person = cost / n_c,
           qaly_per_person = qaly / n_c),
  icer_table_ltbi_rw %>%
    mutate(scenario = "realworld_55pct",
           cost_per_person = cost / n_c,
           qaly_per_person = qaly / n_c)
)
write.csv(icer_table_ltbi_combined, "output/csv/ltbi_completion_sensitivity.csv", row.names = FALSE)
cat("Saved: output/csv/ltbi_completion_sensitivity.csv\n\n")

# Frontier comparison plot: base -> real-world arrows
frontier_names <- icer_table %>%
  filter(dominance %in% c("ref", "non-dominated")) %>%
  arrange(cost) %>%
  pull(strategy)

frontier_ltbi_base <- icer_table %>%
  filter(strategy %in% frontier_names) %>%
  select(strategy, q_base = qaly, c_base = cost) %>%
  mutate(q_base = q_base / n_c, c_base = c_base / n_c)

frontier_ltbi_rw <- icer_table_ltbi_rw %>%
  filter(strategy %in% frontier_names) %>%
  select(strategy, q_rw = qaly, c_rw = cost) %>%
  mutate(q_rw = q_rw / n_c, c_rw = c_rw / n_c)

frontier_ltbi_arrows <- frontier_ltbi_base %>%
  left_join(frontier_ltbi_rw, by = "strategy") %>%
  mutate(strategy = factor(strategy, levels = frontier_names))

cond_label_base <- sprintf("Base (0.253/mo; ~%.0f%% conditional; Surey 2021)", cond_completion_base)
cond_label_rw   <- "Real-world (0.185/mo; 55.5% UKHSA 2025)"

frontier_ltbi_long <- bind_rows(
  frontier_ltbi_arrows %>% transmute(strategy, q = q_base, c = c_base, cond = cond_label_base),
  frontier_ltbi_arrows %>% transmute(strategy, q = q_rw,   c = c_rw,   cond = cond_label_rw)
)

label_ltbi_df <- frontier_ltbi_arrows %>%
  transmute(strategy, q = q_rw, c = c_rw)

p_ltbi_compare <- ggplot() +
  geom_path(data = frontier_ltbi_arrows %>% arrange(q_base),
            aes(x = q_base, y = c_base),
            colour = "#7a5090", linetype = "dashed", linewidth = 0.6, alpha = 0.7) +
  geom_path(data = frontier_ltbi_arrows %>% arrange(q_rw),
            aes(x = q_rw, y = c_rw),
            colour = "#c84b6a", linetype = "solid", linewidth = 0.6, alpha = 0.8) +
  geom_segment(
    data = frontier_ltbi_arrows,
    aes(x = q_base, y = c_base, xend = q_rw, yend = c_rw),
    arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
    colour = "#c84b6a", linewidth = 0.50
  ) +
  geom_point(data = frontier_ltbi_long,
             aes(x = q, y = c, colour = cond, shape = cond),
             size = 4, alpha = 0.95) +
  ggrepel::geom_text_repel(
    data = label_ltbi_df,
    aes(x = q, y = c, label = strategy),
    size = 3.2, colour = "grey20",
    nudge_y = 4, box.padding = 0.3,
    segment.color = "grey60", max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = setNames(c("#7a5090", "#c84b6a"),
                                        c(cond_label_base, cond_label_rw)), name = NULL) +
  scale_shape_manual(values  = setNames(c(16, 17),
                                        c(cond_label_base, cond_label_rw)), name = NULL) +
  theme_minimal(base_size = 13) +
  labs(
    x        = "QALYs per person",
    y        = "Cost per person (\u00a3)",
    title    = NULL,
    subtitle = sprintf("Base (~%.0f%% conditional completion; Surey 2021 3HR) vs real-world (55.5%%; UKHSA 2025)\nArrows show shift in frontier strategies when completion reduced to real-world rate",
                       cond_completion_base)
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 9),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10)
  )

ggsave("output/ltbi_completion_compare.png", p_ltbi_compare,
       width = 11, height = 7, dpi = 300)
cat("Saved: output/ltbi_completion_compare.png\n")

# =============================================================================
# SENSITIVITY ANALYSIS: LTBI PREVALENCE
#
# Base case: 17.8% (Berrocal-Almanza et al. 2022, IGRA-screened UK migrants from
#   high-burden countries; pooled estimate from 11 UK studies)
# Sensitivity: 15.1% (UKHSA TB in England 2025 Prevention report, post-entry
#   migrant screening; more recent national estimate)
#
# Direction of effect: lower LTBI prevalence → fewer LTBI detected →
#   smaller QALY gains from LTBI treatment → frontier ICERs increase
#   (parallel IGRA strategies marginally less attractive, but remain <<£25k)
#
# Init vector adjustment:
#   - Sequential strategies: LatentDiagnosed scaled proportionally (×15100/17800);
#     LatentUndiagnosed = 15100 − LatentDiagnosed; Uninfected = n_c − 1000 − 15100
#   - Parallel IGRA strategies: LatentDiagnosed = 15100 × igra_sensitivity;
#     LatentUndiagnosed = 15100 − LatentDiagnosed; Uninfected = 83900
#   - No_screening: LatentUndiagnosed = 15100; Uninfected = 83900
#   - ActiveUndiagnosed: unchanged at 1% (prev_active; separate parameter)
# =============================================================================
cat("\n================================================================================\n")
cat("    LTBI PREVALENCE SENSITIVITY \u2014 17.8% (base) vs 15.1% (UKHSA 2025)          \n")
cat("================================================================================\n\n")

prev_ltbi_ukhsa    <- as.numeric(config_list["prev_ltbi_ukhsa"])
n_ltbi_ukhsa       <- n_c * prev_ltbi_ukhsa                         # 15,100 per 100k
n_uninfected_ukhsa <- n_c * (1 - prev_ltbi_ukhsa - prev_active)     # 83,900 per 100k
scale_ltbi_prev    <- prev_ltbi_ukhsa / prev_ltbi                    # 15.1 / 17.8 = 0.8483

cat(sprintf("Base:        prev_ltbi = %.1f%% (Berrocal-Almanza 2022) \u2192 %g LTBI per 100k\n",
            prev_ltbi * 100, n_ltbi_global))
cat(sprintf("Sensitivity: prev_ltbi = %.1f%% (UKHSA 2025)             \u2192 %g LTBI per 100k\n\n",
            prev_ltbi_ukhsa * 100, n_ltbi_ukhsa))

strategies_ltbi_prev <- lapply(strategies, function(s) {
  if (is.null(s$excel_name) || is.na(s$excel_name)) {
    # No_screening: no LTBI detected
    s$init["LatentUndiagnosed"] <- n_ltbi_ukhsa
    s$init["Uninfected"]        <- n_uninfected_ukhsa
  } else if (!is.null(s$prog_cost_type)) {
    # Parallel IGRA: detection = n_ltbi_ukhsa × igra_sensitivity
    sens <- if (grepl("tspt", s$excel_name, ignore.case = TRUE)) igra_sens_tspt else igra_sens_qft
    ltbi_det_new <- n_ltbi_ukhsa * sens
    s$init["LatentDiagnosed"]   <- ltbi_det_new
    s$init["LatentUndiagnosed"] <- n_ltbi_ukhsa - ltbi_det_new
    s$init["Uninfected"]        <- n_uninfected_ukhsa
  } else {
    # Sequential strategies: scale LatentDiagnosed proportionally
    ltbi_det_new <- s$init["LatentDiagnosed"] * scale_ltbi_prev
    s$init["LatentDiagnosed"]   <- ltbi_det_new
    s$init["LatentUndiagnosed"] <- n_ltbi_ukhsa - ltbi_det_new
    s$init["Uninfected"]        <- n_uninfected_ukhsa
  }
  s
})

strategy_results_ltbi_prev <- mclapply(strategies_ltbi_prev, run_strategy,
                                        params   = paramsData,
                                        mc.cores = max(1L, detectCores() - 1L))
icer_table_ltbi_prev <- calculate_icer(strategy_results_ltbi_prev)

cat("ICER TABLE \u2014 Sensitivity: LTBI prevalence 15.1% (UKHSA 2025):\n")
cat(sprintf("%-22s %11s %10s %11s %14s %15s  %s\n",
            "Strategy", "Cost/person", "QALYs/pp", "Inc.Cost",
            "ICER(vs ref)", "Sequential ICER", "Dominance"))
cat(strrep("-", 108), "\n")
for (i in 1:nrow(icer_table_ltbi_prev)) {
  icer_str     <- if (is.na(icer_table_ltbi_prev$icer[i])) "Ref"
                  else paste0("\u00a3", format(round(icer_table_ltbi_prev$icer[i]), big.mark = ","))
  seq_icer_str <- if (is.na(icer_table_ltbi_prev$sequential_icer[i])) "-"
                  else paste0("\u00a3", format(round(icer_table_ltbi_prev$sequential_icer[i]),
                                               big.mark = ","))
  cat(sprintf("%-22s %10s %10.4f %10s %14s %15s  %s\n",
              icer_table_ltbi_prev$strategy[i],
              paste0("\u00a3", format(round(icer_table_ltbi_prev$cost_per_person[i]), big.mark = ",")),
              icer_table_ltbi_prev$qaly_per_person[i],
              paste0("\u00a3", format(round(icer_table_ltbi_prev$inc_cost[i] / n_c), big.mark = ",")),
              icer_str,
              seq_icer_str,
              icer_table_ltbi_prev$dominance[i]))
}
cat("\n")

# Frontier comparison: base (17.8%) vs UKHSA 2025 (15.1%)
frontier_names_prev <- icer_table %>%
  filter(dominance %in% c("ref", "non-dominated")) %>%
  arrange(cost) %>%
  pull(strategy)

cat("FRONTIER COMPARISON \u2014 LTBI Prevalence: Base (17.8%) vs UKHSA 2025 (15.1%):\n")
cat(sprintf("  %-32s  %18s  %24s  %12s\n",
    "Strategy", "Base ICER (17.8%)", "UKHSA 2025 ICER (15.1%)", "Delta"))
cat(strrep("-", 96), "\n")
for (st in frontier_names_prev) {
  r_base <- icer_table[icer_table$strategy == st, ]
  r_prev <- icer_table_ltbi_prev[icer_table_ltbi_prev$strategy == st, ]
  icer_b <- if (is.na(r_base$sequential_icer)) "Ref" else
            paste0("\u00a3", format(round(r_base$sequential_icer), big.mark = ","), "/QALY")
  icer_p <- if (is.na(r_prev$sequential_icer)) "Ref" else
            paste0("\u00a3", format(round(r_prev$sequential_icer), big.mark = ","), "/QALY")
  delta  <- if (is.na(r_base$sequential_icer) || is.na(r_prev$sequential_icer)) "\u2014" else
            paste0(ifelse(r_prev$sequential_icer - r_base$sequential_icer >= 0, "+", ""),
                   "\u00a3", format(round(r_prev$sequential_icer - r_base$sequential_icer),
                                    big.mark = ","))
  cat(sprintf("  %-32s  %18s  %24s  %12s\n", st, icer_b, icer_p, delta))
}
cat("\n================================================================================\n\n")

# Save CSV
icer_table_prev_combined <- bind_rows(
  icer_table %>%
    mutate(scenario = "base_17.8pct",
           cost_per_person = cost / n_c,
           qaly_per_person = qaly / n_c),
  icer_table_ltbi_prev %>%
    mutate(scenario = "ukhsa2025_15.1pct",
           cost_per_person = cost / n_c,
           qaly_per_person = qaly / n_c)
)
write.csv(icer_table_prev_combined, "output/csv/ltbi_prevalence_sensitivity.csv", row.names = FALSE)
cat("Saved: output/csv/ltbi_prevalence_sensitivity.csv\n")

# LTBI prevalence sensitivity plot: base (17.8%) vs UKHSA 2025 (15.1%)
frontier_ltbi_prev_base <- icer_table %>%
  filter(strategy %in% frontier_names_prev) %>%
  select(strategy, q_base = qaly, c_base = cost) %>%
  mutate(q_base = q_base / n_c, c_base = c_base / n_c)

frontier_ltbi_prev_new <- icer_table_ltbi_prev %>%
  filter(strategy %in% frontier_names_prev) %>%
  select(strategy, q_new = qaly, c_new = cost) %>%
  mutate(q_new = q_new / n_c, c_new = c_new / n_c)

frontier_prev_arrows <- frontier_ltbi_prev_base %>%
  left_join(frontier_ltbi_prev_new, by = "strategy") %>%
  mutate(strategy = factor(strategy, levels = frontier_names_prev))

lbl_base_prev <- "Base (17.8%; Berrocal-Almanza 2022)"
lbl_new_prev  <- "UKHSA 2025 (15.1%)"

frontier_prev_long <- bind_rows(
  frontier_prev_arrows %>% transmute(strategy, q = q_base, c = c_base, cond = lbl_base_prev),
  frontier_prev_arrows %>% transmute(strategy, q = q_new,  c = c_new,  cond = lbl_new_prev)
)

label_prev_df <- frontier_prev_arrows %>%
  transmute(strategy, q = q_new, c = c_new)

p_ltbi_prev <- ggplot() +
  geom_path(data = frontier_prev_arrows %>% arrange(q_base),
            aes(x = q_base, y = c_base),
            colour = "#7a5090", linetype = "dashed", linewidth = 0.6, alpha = 0.7) +
  geom_path(data = frontier_prev_arrows %>% arrange(q_new),
            aes(x = q_new, y = c_new),
            colour = "#c84b6a", linetype = "solid", linewidth = 0.6, alpha = 0.8) +
  geom_segment(
    data = frontier_prev_arrows,
    aes(x = q_base, y = c_base, xend = q_new, yend = c_new),
    arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
    colour = "#c84b6a", linewidth = 0.50
  ) +
  geom_point(data = frontier_prev_long,
             aes(x = q, y = c, colour = cond, shape = cond),
             size = 4, alpha = 0.95) +
  ggrepel::geom_text_repel(
    data = label_prev_df,
    aes(x = q, y = c, label = strategy),
    size = 3.2, colour = "grey20",
    nudge_y = 4, box.padding = 0.3,
    segment.color = "grey60", max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = setNames(c("#7a5090", "#c84b6a"),
                                        c(lbl_base_prev, lbl_new_prev)), name = NULL) +
  scale_shape_manual(values  = setNames(c(16, 17),
                                        c(lbl_base_prev, lbl_new_prev)), name = NULL) +
  scale_y_continuous(labels = function(x) paste0("\u00a3", x)) +
  theme_minimal(base_size = 13) +
  labs(
    x        = "QALYs per person",
    y        = "Cost per person (\u00a3)",
    subtitle = "Base (17.8%; Berrocal-Almanza 2022) vs UKHSA 2025 (15.1%)\nArrows show shift in frontier strategies at lower LTBI prevalence"
  ) +
  theme(
    plot.subtitle    = element_text(color = "grey40", size = 9),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10)
  )

ggsave("output/ltbi_prevalence_sensitivity.png", p_ltbi_prev,
       width = 11, height = 7, dpi = 300)
cat("Saved: output/ltbi_prevalence_sensitivity.png\n")

# =============================================================================
# ACTIVE TB PREVALENCE SENSITIVITY ANALYSIS
# Base: 1% (1,000/100k; Zenner et al. 2025 ERJ decision tree, high-risk migrants)
# SA1:  0.215% (215/100k; Osei-Yeboah et al. 2025 pooled meta-analysis: 40M+ migrants)
# SA2:  0.44%  (440/100k; upper bound for refugee populations; Osei-Yeboah 2025)
# Method: TP and FN scaled proportionally (fixed detection rate assumption).
# Rationale: entry prevalence does not alter test accuracy; it scales the absolute
#   number of cases detected and missed proportionally.
# =============================================================================
cat("\n================================================================================\n")
cat("         ACTIVE TB PREVALENCE \u2014 ENTRY PREVALENCE SENSITIVITY                    \n")
cat("================================================================================\n\n")
cat(sprintf("Base case:   prev_active = %.1f%% (%d/100k)  [Zenner et al. 2025 ERJ]\n",
            prev_active * 100, round(n_c * prev_active)))
cat("Sensitivity: 215/100k  [Osei-Yeboah 2025 pooled meta-analysis: 40M+ migrants]\n")
cat("             440/100k  [Osei-Yeboah 2025: high-risk subgroup estimate]\n\n")

active_tb_sa_scenarios <- list(
  list(label = "base_1pct",      prev = prev_active,                                      desc = "1%   (1,000/100k; base case)"),
  list(label = "pooled_0215pct", prev = as.numeric(config_list["prev_active_low"]),        desc = "0.215% (215/100k; Osei-Yeboah 2025 pooled)"),
  list(label = "highrisk_044pct", prev = as.numeric(config_list["prev_active_high"]),      desc = "0.44%  (440/100k; high-risk subgroup; Osei-Yeboah et al. 2025)")
)

active_tb_sa_results <- list()

for (sc in active_tb_sa_scenarios) {
  scale_f <- sc$prev / prev_active
  strats_atb <- lapply(strategies, function(s) {
    s2 <- s
    if (!is.null(s2$base_tp) && !is.na(s2$base_tp)) {
      new_tp <- s2$base_tp * scale_f
      new_fn <- n_c * sc$prev - new_tp
      s2$init["ActiveDiagnosed"]   <- max(0, new_tp)
      s2$init["ActiveUndiagnosed"] <- max(0, new_fn)
    } else {
      # No_screening
      s2$init["ActiveDiagnosed"]   <- 0
      s2$init["ActiveUndiagnosed"] <- n_c * sc$prev
    }
    s2
  })
  res_atb <- mclapply(strats_atb, run_strategy,
                      params   = paramsData,
                      mc.cores = max(1L, detectCores() - 1L))
  icer_atb <- calculate_icer(res_atb)
  active_tb_sa_results[[sc$label]] <- list(icer = icer_atb, desc = sc$desc, prev = sc$prev)
}

# Print frontier comparison
cat("FRONTIER COMPARISON \u2014 Active TB Entry Prevalence:\n")
frontier_strategies <- active_tb_sa_results[["base_1pct"]]$icer %>%
  filter(dominance == "non-dominated") %>%
  pull(strategy)
frontier_strategies <- c("No screening", frontier_strategies[frontier_strategies != "No screening"])

header <- sprintf("  %-40s  %12s  %12s  %12s",
                  "Strategy",
                  "Base 1%",
                  "Pooled 0.215%",
                  "Refugee 0.44%")
cat(header, "\n")
cat(strrep("-", nchar(header)), "\n")

for (strat in frontier_strategies) {
  icers <- sapply(active_tb_sa_results, function(sc) {
    r <- sc$icer %>% filter(strategy == strat)
    if (nrow(r) == 0) return("N/A")
    si <- r$sequential_icer
    if (is.na(si) || si == Inf) return("dominated")
    sprintf("\u00a3%s/QALY", format(round(si), big.mark = ","))
  })
  cat(sprintf("  %-40s  %12s  %12s  %12s\n", strat, icers[1], icers[2], icers[3]))
}

# Save CSV
atb_sa_all <- bind_rows(lapply(names(active_tb_sa_results), function(nm) {
  active_tb_sa_results[[nm]]$icer %>%
    mutate(scenario       = nm,
           active_tb_prev = active_tb_sa_results[[nm]]$prev)
}))
write.csv(atb_sa_all, "output/csv/active_tb_prev_sensitivity.csv", row.names = FALSE)
cat("\nSaved: output/csv/active_tb_prev_sensitivity.csv\n")

# Simple comparison plot: ICER for frontier strategies across prevalence scenarios
atb_plot_df <- bind_rows(lapply(names(active_tb_sa_results), function(nm) {
  sc <- active_tb_sa_results[[nm]]
  sc$icer %>%
    filter(strategy %in% frontier_strategies, strategy != "No screening") %>%
    mutate(
      scenario = nm,
      prev_pct = sc$prev * 100,
      label    = sc$desc
    )
})) %>%
  filter(!is.na(sequential_icer) & sequential_icer != Inf & sequential_icer > 0)

if (nrow(atb_plot_df) > 0) {
  # x-axis order: ascending prevalence (0.215 < 0.44 < 1)
  atb_plot_df$prev_label <- factor(
    sprintf("%.3g%%", atb_plot_df$prev_pct),
    levels = c("0.215%", "0.44%", "1%")
  )

  p_atb_sa <- ggplot(atb_plot_df,
                     aes(x = prev_label, y = sequential_icer, fill = strategy)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
    # Annotate missing Cough+CXR bar at 0.215% (extendedly dominated — parallel
    # IGRA strategies produce lower incremental cost per QALY at low active TB burden)
    annotate("text", x = 1, y = 2000,
             label = "Cough+CXR\next. dominated",
             colour = "#9870b8", size = 3, hjust = 0.5, fontface = "italic") +
    geom_hline(yintercept = 25000, linetype = "dashed", colour = "red", linewidth = 0.8) +
    annotate("text", x = 2.6, y = 26500, label = "\u00a325,000/QALY (NICE threshold)",
             colour = "red", size = 3.5, hjust = 1) +
    scale_y_continuous(labels = scales::label_comma(prefix = "\u00a3"),
                       expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = c("Cough+CXR (TB sx)"            = "#9870b8",
                                  "Parallel Cough+QFT (Ultra)"   = "#c84b6a",
                                  "Parallel Sx+QFT (Ultra)"      = "#7a1a4a")) +
    labs(
      title    = NULL,
      subtitle = paste0("Sequential ICERs for frontier strategies | Cough+CXR ext. dominated at 0.215% prevalence\n",
                        "(insufficient active TB cases for sequential screening to be CE vs parallel IGRA)"),
      x        = "Active TB entry prevalence (%)",
      y        = "Sequential ICER (\u00a3/QALY)",
      fill     = "Strategy"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey40"),
          legend.position = "bottom")
  ggsave("output/active_tb_prev_sensitivity.png", p_atb_sa,
         width = 10, height = 6, dpi = 300)
  cat("Saved: output/active_tb_prev_sensitivity.png\n")
}

cat("\n================================================================================\n")

# =============================================================================
# SCENARIO ANALYSIS: COST OF UNMANAGED TB
#
# The universal simplifying assumption in published UK TB CEA models is that
# undiagnosed and not-treated active TB states carry zero NHS cost per cycle
# (Jit 2011 BMJ PMC3273731; Dale 2022 AJE PMID34017976; Green 2025 ERJ PMC12183743).
# This is a known limitation: symptomatic TB patients make repeated NHS contact
# before diagnosis (~6-8 week median delay; UKHSA 2024), and non-adherent patients
# incur Enhanced Case Management costs under NICE NG33.
#
# This scenario is a novel sensitivity check: we assign plausible per-cycle NHS
# contact costs to both unmanaged active TB states, built from UK unit costs.
#
#   Cost_ActiveUndiagnosed_scenario = 150/mo
#     Derived from UK pre-diagnosis contact data:
#     — Loutet et al. 2018 (PHE, England, n=22,422; PMID 29923481): median
#       TOTAL delay 2.8 months (patient 1.3mo + healthcare 14 days; pulmonary TB)
#     — Mawer et al. 2007 (UK BJGP; PMID 17263928): 2-4 GP visits before dx
#       = 0.7-1.4/mo over 2.8 months x 49 (PSSRU 2022) = 34-69/mo
#     — Schwartzman et al. 2002 (Canada; PMID 11936743): 47% have >=1 ED visit
#       pre-diagnosis, mean 2.2 visits over 2.8 months = 0.37/mo x 180
#       (NHS NCC 2022/23) = 67/mo
#     — First-principles range: 100-135/mo; 150 adopted as conservative midpoint
#     — Methodological precedent: Miners et al. 2017 Lancet HIV (UK; PMC5614770)
#       assigns non-zero pre-diagnosis costs in a UK infectious disease Markov model
#     — LOW CONFIDENCE: no UK study directly quantifies pre-diagnosis TB NHS
#       contact costs; this is a model assumption built from unit costs
#
#   Cost_ActiveNotreated_scenario = 275/mo  [LOW CONFIDENCE — model estimate]
#     — Symptom-driven contacts same as ActiveUndiagnosed (~150/mo): known TB
#       patient with progressing untreated disease attends GP/A&E at same or
#       higher rate (Loutet 2018; Mawer 2007; Schwartzman 2002)
#     — Plus ECM outreach (~125/mo): NICE NG33 Enhanced Case Management;
#       Hayward/Holden 2019 Lancet (PMID 30799062): DOT 44/visit; refusing
#       patients receive Level 1-2 ECM = weekly/fortnightly home visits;
#       Hanif 2017 BMC PH (PMID 29141600): Level 1-2 = weekly/fortnightly
#       = ~3 visits/month x 44 = 132/mo; rounded to 125/mo
#     — No published data on contact frequency specifically for refusing patients;
#       UKHSA records refusal rates (~1-2%) but not contact frequencies
#
# Ref: Loutet 2018 PMID29923481; Mawer 2007 PMID17263928; Schwartzman 2002
#      PMID11936743; NHS NCC 2022/23; PSSRU 2022; NICE NG33; Miners 2017
#      PMC5614770; Hayward 2019 PMID30799062; Hanif 2017 PMID29141600.
# =============================================================================
cat("\n")
cat("================================================================================\n")
cat("         SCENARIO ANALYSIS: COST OF UNMANAGED TB                               \n")
cat("================================================================================\n\n")

inaction_ActiveUndiagnosed <- as.numeric(config_list["Cost_ActiveUndiagnosed_scenario"])
inaction_ActiveNotreated   <- as.numeric(config_list["Cost_ActiveNotreated_scenario"])

cat(sprintf("Base case:  Cost_ActiveUndiagnosed = £0/mo  |  Cost_ActiveNotreated = £0/mo\n"))
cat(sprintf("Scenario:   Cost_ActiveUndiagnosed = £%g/mo  |  Cost_ActiveNotreated = £%g/mo\n",
            inaction_ActiveUndiagnosed, inaction_ActiveNotreated))
cat(sprintf("(Novel sensitivity — zero cost is universal published convention; see comment block above)\n\n"))

params_inaction <- paramsData
params_inaction[["Cost_ActiveUndiagnosed"]] <- inaction_ActiveUndiagnosed
params_inaction[["Cost_ActiveNotreated"]]   <- inaction_ActiveNotreated

strategy_results_inaction <- mclapply(strategies, run_strategy,
                                      params     = params_inaction,
                                      mc.cores   = max(1L, detectCores() - 1L))
icer_table_inaction        <- calculate_icer(strategy_results_inaction)

cat("ICER TABLE — Scenario: introducing costs for unmanaged active TB (base case: \u00a30):\n")
cat(sprintf("%-22s %11s %10s %11s %14s %15s  %s\n",
            "Strategy", "Cost/person", "QALYs/pp", "Inc.Cost",
            "ICER(vs ref)", "Sequential ICER", "Dominance"))
cat(strrep("-", 108), "\n")
for (i in 1:nrow(icer_table_inaction)) {
  icer_str     <- if (is.na(icer_table_inaction$icer[i])) "Ref"
                  else paste0("£", format(round(icer_table_inaction$icer[i]), big.mark = ","))
  seq_icer_str <- if (is.na(icer_table_inaction$sequential_icer[i])) "-"
                  else paste0("£", format(round(icer_table_inaction$sequential_icer[i]), big.mark = ","))
  cat(sprintf("%-22s %10s %10.4f %10s %14s %15s  %s\n",
              icer_table_inaction$strategy[i],
              paste0("£", format(round(icer_table_inaction$cost_per_person[i]), big.mark = ",")),
              icer_table_inaction$qaly_per_person[i],
              paste0("£", format(round(icer_table_inaction$inc_cost[i] / n_c), big.mark = ",")),
              icer_str,
              seq_icer_str,
              icer_table_inaction$dominance[i]))
}
cat("\n================================================================================\n\n")

write.csv(icer_table_inaction %>%
  mutate(cost_per_person = cost / n_c, qaly_per_person = qaly / n_c),
  "output/csv/icer_table_inaction.csv", row.names = FALSE)
cat("Saved: output/csv/icer_table_inaction.csv\n")

# Frontier-only arrow plot: base case -> cost of unmanaged TB scenario
# Shows only the 4 frontier strategies with arrows indicating direction of shift.
# Cleaner than full-35 scatter; makes the cost increase for No Screening visually obvious.

# Derive frontier names dynamically from ICER table (non-dominated + ref)
frontier_names <- icer_table %>%
  filter(dominance %in% c("ref", "non-dominated")) %>%
  arrange(cost) %>%
  pull(strategy)

frontier_base <- icer_table %>%
  filter(strategy %in% frontier_names) %>%
  select(strategy, q_base = qaly, c_base = cost) %>%
  mutate(q_base = q_base / n_c, c_base = c_base / n_c)

frontier_scen <- icer_table_inaction %>%
  filter(strategy %in% frontier_names) %>%
  select(strategy, q_scen = qaly, c_scen = cost) %>%
  mutate(q_scen = q_scen / n_c, c_scen = c_scen / n_c)

frontier_arrows <- left_join(frontier_base, frontier_scen, by = "strategy") %>%
  mutate(strategy = factor(strategy, levels = frontier_names))

# Long form for points + legend
frontier_long <- bind_rows(
  frontier_arrows %>% transmute(strategy, q = q_base, c = c_base,
                                 cond = "Base case (\u00a30 for unmanaged TB)"),
  frontier_arrows %>% transmute(strategy, q = q_scen, c = c_scen,
                                 cond = "Scenario (\u00a3150\u2013275/mo for unmanaged TB)")
)

# Label positions: use scenario (higher) point, nudge up
label_df <- frontier_arrows %>%
  transmute(strategy, q = q_scen, c = c_scen)

p_inaction_compare <- ggplot() +
  # Frontier line — base case (dashed)
  geom_path(data = frontier_arrows %>% arrange(q_base),
            aes(x = q_base, y = c_base),
            colour = "#7a5090", linetype = "dashed", linewidth = 0.6, alpha = 0.7) +
  # Frontier line — scenario (solid)
  geom_path(data = frontier_arrows %>% arrange(q_scen),
            aes(x = q_scen, y = c_scen),
            colour = "#c84b6a", linetype = "solid", linewidth = 0.6, alpha = 0.7) +
  # Arrows: base -> scenario for each frontier strategy
  geom_segment(
    data = frontier_arrows,
    aes(x = q_base, y = c_base, xend = q_scen, yend = c_scen),
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    colour = "grey40", linewidth = 0.5
  ) +
  # Points
  geom_point(data = frontier_long,
             aes(x = q, y = c, colour = cond, shape = cond),
             size = 4, alpha = 0.95) +
  # Strategy labels (above scenario point)
  ggrepel::geom_text_repel(
    data = label_df,
    aes(x = q, y = c, label = strategy),
    size = 3.2, colour = "grey20",
    nudge_y = 4, box.padding = 0.3,
    segment.color = "grey60", max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("Base case (\u00a30 for unmanaged TB)" = "#7a5090",
               "Scenario (\u00a3150\u2013275/mo for unmanaged TB)" = "#c84b6a"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("Base case (\u00a30 for unmanaged TB)" = 16,
               "Scenario (\u00a3150\u2013275/mo for unmanaged TB)" = 17),
    name = NULL
  ) +
  theme_minimal(base_size = 13) +
  labs(
    x        = "QALYs per person",
    y        = "Cost per person (\u00a3)",
    title    = NULL,
    subtitle = "Arrows show direction of shift from base case (circles) to scenario (triangles) | Frontier strategies only"
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11)
  )

ggsave("output/scenario_inaction_compare.png", p_inaction_compare,
       width = 10, height = 7, dpi = 300)
cat("Saved: output/scenario_inaction_compare.png\n")

cat("\n================================================================================\n")

# =============================================================================
# SCENARIO 2: TRANSMISSION PREVENTION
#
# Base case assumes no onward transmission from undiagnosed active TB cases
# (standard convention in UK TB CEA models; Jit 2011, Dale 2022, Green 2025).
# This scenario quantifies the value of averted secondary transmissions.
#
# Method (Brooks-Pollock et al. 2020, PLOS Comput Biol, PMID 32218567):
#   1. Estimate new active TB cases per strategy over the model horizon:
#      new_atb = p_reactivation × cumulative person-months in all reactivating
#      latent states (LatentUndiagnosed, LatentDiagnosed, LatentNotreated,
#      LatentLtfu at 0.000912/mo; LatentDiscontinued at 0.0008/mo).
#   2. Active TB cases prevented vs No Screening = no_screening_atb - strategy_atb.
#   3. Secondary cases avoided = cases_prevented × beta_tx
#      (base: beta = 0.205, conservative base = 5 contacts/case/yr × 4.1%
#       transmission probability per contact; anchored to Brooks-Pollock et al.
#       2020 UK estimate (β=0.41 at 10 contacts);
#       DSA: 5 (base)/10 (mid)/15 (high) contacts/yr; PSA: contacts ~ Uniform(5,15)).
#   4. Cost saving = secondary cases avoided × cost per secondary TB case.
#      (Green et al. 2025 ERJ Open Res: £6,055 base; DSA range £3,028–£12,110.)
#   5. Subtract savings from strategy costs → recalculate ICERs.
#
# This approach is consistent with Pareek et al. 2011 (Lancet ID) and
# Green et al. 2025 (ERJ Open Res) which apply a fixed secondary-case
# multiplier to active TB cases prevented.
#
# Outputs: icer_table_transmission.csv, scenario_transmission_compare.png
# =============================================================================
cat("\n")
cat("================================================================================\n")
cat("         SCENARIO ANALYSIS: TRANSMISSION PREVENTION                            \n")
cat("================================================================================\n\n")

p_tx_per_contact      <- as.numeric(config_list["p_tx_per_contact"])
n_contacts_base       <- as.numeric(config_list["n_contacts_base"])
n_contacts_dsa_mid    <- as.numeric(config_list["n_contacts_dsa_mid"])
n_contacts_dsa_high   <- as.numeric(config_list["n_contacts_dsa_high"])
cost_per_secondary_tb <- as.numeric(config_list["cost_secondary_tb"])
beta_tx_base          <- n_contacts_base * p_tx_per_contact
p_react_standard      <- 0.000912  # monthly reactivation rate phase 1 (Berrocal-Almanza 2022, config row 41)
p_react_phase2_tx     <- as.numeric(config_list["p_react_phase2"])  # 0.000333, phase 2 rate
p_react_disc          <- 0.0008    # LatentDiscontinued reactivation (config row 74)
# Scale LatentDiscontinued rate for phase 2 proportionally (same phase switch logic as main model)
p_react_disc_phase2   <- p_react_disc * (p_react_phase2_tx / p_react_standard)

cat(sprintf("Beta (secondary cases per active TB case, UK): %.2f (Brooks-Pollock 2020, 95%% CrI 0.30-0.60)\n", beta_tx_base))
cat(sprintf("Contacts parameterisation: %d contacts/case/yr × %.3f transmission probability = beta %.2f\n",
            n_contacts_base, p_tx_per_contact, beta_tx_base))
cat(sprintf("Cost per secondary case:  £%s (Green et al. 2025 ERJ Open Res, PMC12183743; NHS direct costs, 2024 GBP)\n",
            format(cost_per_secondary_tb, big.mark = ",")))
cat(sprintf("(DSA: contacts = %g / %g / %g per untreated case/yr; PSA: Uniform(%g, %g))\n\n",
            n_contacts_base, n_contacts_dsa_mid, n_contacts_dsa_high,
            n_contacts_base, n_contacts_dsa_high))

# Step 1: Estimate new active TB cases per strategy over the 55-year horizon.
# New active TB cases arise from reactivation out of all reactivating latent states.
# Approximation: new_cases_from_state = p_reactivation × cumulative person-months in state
# (standard Markov cohort incidence approximation; cycle length = 1 month).
# Two-phase rates applied: phase 1 (cycles 1–60) = 0.000912/mo;
#                          phase 2 (cycles 61–660) = 0.000333/mo (Horsburgh 2004 NEJM).
latent_states_standard <- c("LatentUndiagnosed", "LatentDiagnosed",
                             "LatentNotreated",  "LatentLtfu")

new_atb_cases <- sapply(strategy_results, function(x) {
  sm      <- x$state_membership
  n_cy    <- nrow(sm)
  cut     <- min(60L, n_cy)         # phase 1: cycles 1-60
  # Phase 1 person-months
  std_p1  <- sum(sm[seq_len(cut),          latent_states_standard])
  disc_p1 <- sum(sm[seq_len(cut),          "LatentDiscontinued"])
  # Phase 2 person-months (cycles 61+)
  std_p2  <- if (n_cy > cut) sum(sm[(cut + 1L):n_cy, latent_states_standard]) else 0
  disc_p2 <- if (n_cy > cut) sum(sm[(cut + 1L):n_cy, "LatentDiscontinued"])   else 0

  std_p1  * p_react_standard   + disc_p1 * p_react_disc +
  std_p2  * p_react_phase2_tx  + disc_p2 * p_react_disc_phase2
})
names(new_atb_cases) <- sapply(strategy_results, function(x) x$strategy_name)

# Step 2: Active TB cases prevented vs No Screening
atb_no_screening  <- new_atb_cases["No screening"]
atb_cases_prevent <- atb_no_screening - new_atb_cases   # positive = screening prevents

# Step 3: Secondary cases avoided = active TB cases prevented × beta
secondary_cases_avoided <- atb_cases_prevent * beta_tx_base

# Step 4: Total cost savings for the 100k cohort
transmission_savings <- secondary_cases_avoided * cost_per_secondary_tb

cat("New active TB cases (from reactivation) and transmission savings by strategy:\n")
cat(sprintf("  %-28s %14s %16s %18s %16s\n",
            "Strategy", "New ATB (100k)", "ATB prevented", "Sec. cases avoid.", "Savings (£, 100k)"))
cat(strrep("-", 96), "\n")
for (nm in names(new_atb_cases)) {
  cat(sprintf("  %-28s %14.0f %16.0f %18.1f %16s\n",
              nm,
              new_atb_cases[nm],
              atb_cases_prevent[nm],
              secondary_cases_avoided[nm],
              format(round(transmission_savings[nm]), big.mark = ",")))
}
cat("\n")

# Save transmission prevention counts to CSV
tx_counts_df <- data.frame(
  strategy              = names(new_atb_cases),
  new_atb_cases         = as.numeric(new_atb_cases),
  atb_cases_prevented   = as.numeric(atb_cases_prevent),
  secondary_cases_avoided = as.numeric(secondary_cases_avoided),
  transmission_saving_total = as.numeric(transmission_savings),
  transmission_saving_pp    = as.numeric(transmission_savings) / n_c,
  stringsAsFactors = FALSE
)
write.csv(tx_counts_df, "output/csv/transmission_prevention_counts.csv", row.names = FALSE)
cat("Saved: output/csv/transmission_prevention_counts.csv\n")

# Step 5: Build modified strategy_results with adjusted total_cost, then recalculate ICERs
strategy_results_transmission <- lapply(strategy_results, function(x) {
  x_mod <- x
  savings <- transmission_savings[x$strategy_name]
  if (!is.na(savings)) {
    x_mod$total_cost <- x$total_cost - savings
  }
  x_mod
})

icer_table_transmission <- calculate_icer(strategy_results_transmission)

cat(sprintf("ICER TABLE — Scenario 2: Secondary Cases Prevention (beta = %.2f, Brooks-Pollock 2020):\n", beta_tx_base))
cat(sprintf("%-22s %11s %10s %11s %14s %15s  %s\n",
            "Strategy", "Cost/person", "QALYs/pp", "Inc.Cost",
            "ICER(vs ref)", "Sequential ICER", "Dominance"))
cat(strrep("-", 108), "\n")
for (i in 1:nrow(icer_table_transmission)) {
  icer_str     <- if (is.na(icer_table_transmission$icer[i])) "Ref"
                  else paste0("\u00a3", format(round(icer_table_transmission$icer[i]), big.mark = ","))
  seq_icer_str <- if (is.na(icer_table_transmission$sequential_icer[i])) "-"
                  else paste0("\u00a3", format(round(icer_table_transmission$sequential_icer[i]), big.mark = ","))
  cat(sprintf("%-22s %10s %10.4f %10s %14s %15s  %s\n",
              icer_table_transmission$strategy[i],
              paste0("\u00a3", format(round(icer_table_transmission$cost_per_person[i]), big.mark = ",")),
              icer_table_transmission$qaly_per_person[i],
              paste0("\u00a3", format(round(icer_table_transmission$inc_cost[i] / n_c), big.mark = ",")),
              icer_str,
              seq_icer_str,
              icer_table_transmission$dominance[i]))
}
cat("\n================================================================================\n\n")

write.csv(icer_table_transmission %>%
  mutate(cost_per_person = cost / n_c, qaly_per_person = qaly / n_c),
  "output/csv/icer_table_transmission.csv", row.names = FALSE)
cat("Saved: output/csv/icer_table_transmission.csv\n")

# Frontier-only arrow plot: base case → transmission prevention scenario
# Shows frontier strategies with horizontal arrows (costs shift left = cheaper).
# QALYs unchanged; only costs change in this scenario.

# Union of base-case + transmission frontier strategies:
# Frontier composition is unchanged between base case and Scenario 2
# (same 4 strategies: No screening, Cough+CXR, Symptom screen+CXR, Parallel Sx+QFT).
# Union taken dynamically; no strategies added or dropped.
tx_frontier_names <- union(
  frontier_names,
  icer_table_transmission %>%
    filter(dominance %in% c("ref", "non-dominated")) %>%
    arrange(cost) %>% pull(strategy)
)

frontier_base_tx <- icer_table %>%
  filter(strategy %in% tx_frontier_names) %>%
  select(strategy, q_base = qaly, c_base = cost) %>%
  mutate(q_base = q_base / n_c, c_base = c_base / n_c)

frontier_scen_tx <- icer_table_transmission %>%
  filter(strategy %in% tx_frontier_names) %>%
  select(strategy, q_scen = qaly, c_scen = cost) %>%
  mutate(q_scen = q_scen / n_c, c_scen = c_scen / n_c)

frontier_arrows_tx <- left_join(frontier_base_tx, frontier_scen_tx, by = "strategy") %>%
  mutate(strategy = factor(strategy, levels = tx_frontier_names))

lbl_base_tx <- "Base case (no secondary cases)"
lbl_scen_tx <- sprintf("Scenario 2: secondary cases included (\u03b2 = %.2f; \u00a3%s/secondary case)",
                        beta_tx_base,
                        format(cost_per_secondary_tb, big.mark = ","))

frontier_long_tx <- bind_rows(
  frontier_arrows_tx %>% transmute(strategy, q = q_base, c = c_base,
                                   cond = lbl_base_tx),
  frontier_arrows_tx %>% transmute(strategy, q = q_scen, c = c_scen,
                                   cond = lbl_scen_tx)
)

label_df_tx <- frontier_arrows_tx %>%
  transmute(strategy, q = q_base, c = c_base,
            is_igra = grepl("^Parallel", strategy, ignore.case = FALSE))

tx_colours <- c("#7a5090", "#9b5fc0")
names(tx_colours) <- c(lbl_base_tx, lbl_scen_tx)
tx_shapes <- c(16L, 17L)
names(tx_shapes) <- c(lbl_base_tx, lbl_scen_tx)

p_transmission_compare <- ggplot() +
  geom_path(data = frontier_arrows_tx %>% arrange(q_base),
            aes(x = q_base, y = c_base),
            colour = "#7a5090", linetype = "dashed", linewidth = 0.6, alpha = 0.7) +
  geom_path(data = frontier_arrows_tx %>% arrange(q_scen),
            aes(x = q_scen, y = c_scen),
            colour = "#9b5fc0", linetype = "solid", linewidth = 0.6, alpha = 0.7) +
  geom_segment(
    data = frontier_arrows_tx,
    aes(x = q_base, y = c_base, xend = q_scen, yend = c_scen),
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    colour = "grey40", linewidth = 0.5
  ) +
  geom_point(data = frontier_long_tx,
             aes(x = q, y = c, colour = cond, shape = cond),
             size = 4, alpha = 0.95) +
  # IGRA strategy labels in purple-bold; non-IGRA in grey
  ggrepel::geom_text_repel(
    data = label_df_tx %>% filter(!is_igra),
    aes(x = q, y = c, label = strategy),
    size = 3.2, colour = "grey30",
    nudge_y = -8, box.padding = 0.3,
    segment.color = "grey60", max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = label_df_tx %>% filter(is_igra),
    aes(x = q, y = c, label = strategy),
    size = 3.2, colour = "#6a2a90", fontface = "bold",
    nudge_x = 0.004, nudge_y = c(18, 6),
    box.padding = 0.4, point.padding = 0.3,
    segment.color = "#9a6ab0", max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = tx_colours, name = NULL) +
  scale_shape_manual(values  = tx_shapes,  name = NULL) +
  theme_minimal(base_size = 13) +
  labs(
    x        = "QALYs per person",
    y        = "Cost per person (\u00a3)",
    title    = NULL,
    subtitle = sprintf(
      "Arrows show cost reduction when secondary cases are included | Frontier strategies only\n\u03b2 = %.2f secondary cases per active TB case | Cost per secondary case: \u00a3%s\n\u25b2 Parallel IGRA strategies (bold purple) benefit most — more LTBI detected \u2192 more secondary cases prevented",
      beta_tx_base, format(cost_per_secondary_tb, big.mark = ","))
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10)
  )

ggsave("output/scenario_transmission_compare.png", p_transmission_compare,
       width = 10, height = 7, dpi = 300)
cat("Saved: output/scenario_transmission_compare.png\n")

# ---- 2-panel scenario comparison (A = Sc.1, B = Sc.2) -----------------------
# Shared axis limits so panels are directly comparable
panel_xlim <- c(21.455, 21.515)
panel_ylim <- c(270, 480)

p_panel_sc1 <- p_inaction_compare +
  coord_cartesian(xlim = panel_xlim, ylim = panel_ylim) +
  labs(
    title    = NULL,
    subtitle = "Arrows show shift from base case (\u25cf) to scenario (\u25b2) | Frontier strategies only",
    x = "QALYs per person", y = "Cost per person (\u00a3)"
  ) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "grey40", size = 9),
        legend.position = "bottom", legend.text = element_text(size = 9))

# Rebuild IGRA labels in panel B with better separation
p_panel_sc2 <- p_transmission_compare +
  coord_cartesian(xlim = panel_xlim, ylim = panel_ylim) +
  labs(
    title    = NULL,
    subtitle = "Arrows show cost reduction when secondary cases included (\u03b2 = 0.41; \u00a3 6,055/secondary case)",
    x = "QALYs per person", y = ""
  ) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "grey40", size = 9),
        legend.position = "bottom", legend.text = element_text(size = 9))

p_scenario_panel <- (p_panel_sc1 | p_panel_sc2) +
  patchwork::plot_annotation(
    title   = NULL,
    caption = "Sc. 1: unmanaged active TB assigned costs (£150/mo undiagnosed; £275/mo untreated). Sc. 2: secondary cases included (\u03b2 = 0.41; £6,055/secondary case).",
    theme   = theme(plot.title   = element_text(face = "bold", size = 14),
                    plot.caption = element_text(color = "grey40", size = 9))
  )

ggsave("output/scenario_comparison_panel.png", p_scenario_panel,
       width = 16, height = 7, dpi = 300)
cat("Saved: output/scenario_comparison_panel.png\n")

# ---- Dumbbell chart: ICER vs No Screening — base case vs Scenario 2 ----------
dumbbell_strats <- c("Cough+CXR (TB sx)", "Symptom screen+CXR", "Parallel Sx+QFT (Ultra)")

db_base <- icer_table %>%
  filter(strategy %in% dumbbell_strats) %>%
  transmute(strategy, icer_base = icer,
            dominated_base = dominance == "extendedly dominated")

db_sc2 <- icer_table_transmission %>%
  filter(strategy %in% dumbbell_strats) %>%
  transmute(strategy, icer_sc2 = icer)

db_df <- left_join(db_base, db_sc2, by = "strategy") %>%
  mutate(
    is_igra  = grepl("^Parallel", strategy),
    strategy = factor(strategy, levels = rev(dumbbell_strats))
  )

p_dumbbell <- ggplot(db_df, aes(y = strategy)) +
  # connecting line
  geom_segment(aes(x = icer_base, xend = icer_sc2,
                   yend = strategy,
                   colour = is_igra),
               linewidth = 1.8, alpha = 0.6) +
  # base case dot
  geom_point(aes(x = icer_base), shape = 21,
             fill = "white", colour = "grey40", size = 4, stroke = 1.2) +
  # Sc.2 dot
  geom_point(aes(x = icer_sc2, fill = is_igra),
             shape = 21, colour = "white", size = 5) +
  # label: base value (below dot to avoid overlap with Sc.2 label)
  geom_text(aes(x = icer_base,
                label = paste0("\u00a3", format(round(icer_base), big.mark = ","))),
            vjust = 2.2, size = 3, colour = "grey50") +
  # label: Sc.2 value (above dot)
  geom_text(aes(x = icer_sc2,
                label = paste0("\u00a3", format(round(icer_sc2), big.mark = ",")),
                colour = is_igra),
            vjust = -1.2, size = 3.2, fontface = "bold") +
  # ext. dominated annotation
  geom_text(data = db_df %>% filter(dominated_base),
            aes(x = icer_base, label = "(ext. dominated\nin base case)"),
            vjust = 2, hjust = 0.5, size = 2.8, colour = "grey50", fontface = "italic") +
  scale_colour_manual(values = c("FALSE" = "#c84b6a", "TRUE" = "#9b5fc0"), guide = "none") +
  scale_fill_manual(values   = c("FALSE" = "#c84b6a", "TRUE" = "#9b5fc0"), guide = "none") +
  scale_x_continuous(labels = function(x) paste0("\u00a3", format(x, big.mark = ",")),
                     limits = c(900, 9000)) +
  theme_minimal(base_size = 13) +
  labs(
    x        = "ICER vs No Screening (\u00a3/QALY)",
    y        = "",
    title    = NULL,
    subtitle = "Open circle = base case | Filled circle = Scenario 2 (5 contacts/case; conservative)\nFrontier composition unchanged; Parallel Sx+QFT ICER falls ~\u00a3785/QALY from secondary case savings"
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(colour = "grey40", size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y      = element_text(size = 12)
  )

ggsave("output/sc2_dumbbell.png", p_dumbbell, width = 10, height = 5, dpi = 300)
cat("Saved: output/sc2_dumbbell.png\n")

cat("\n================================================================================\n")

# =============================================================================
# TRANSMISSION DSA: contacts/case/yr (5/10/15) × cost (£6,055/£12,110)
#
# Contacts DSA range: 5 (conservative) / 10 (base) / 15 (high) per untreated
# active TB case per year; consistent with WHO/CDC contact investigation literature.
# Implied beta: contacts × p_tx_per_contact (0.041), anchored so that
# 10 contacts × 0.041 = beta 0.41 (Brooks-Pollock 2020 UK estimate).
# cost DSA: £6,055 base vs £12,110 upper (Green 2025 ERJ Open Res DSA range).
# Active TB cases per strategy held fixed from base case run above.
# =============================================================================

tx_dsa_contacts <- c(n_contacts_base, n_contacts_dsa_mid, n_contacts_dsa_high)
tx_dsa_results <- lapply(tx_dsa_contacts, function(nc) {
  beta_eff <- nc * p_tx_per_contact
  savings  <- atb_cases_prevent * beta_eff * cost_per_secondary_tb
  sr_mod   <- lapply(strategy_results, function(x) {
    x$total_cost <- x$total_cost - savings[x$strategy_name]
    x
  })
  tbl <- calculate_icer(sr_mod)
  tbl$n_contacts <- nc
  tbl$beta_eff   <- beta_eff
  tbl
})
tx_dsa_df <- bind_rows(tx_dsa_results)

# Save DSA transmission table
write.csv(tx_dsa_df %>% mutate(cost_per_person = cost / n_c, qaly_per_person = qaly / n_c),
          "output/csv/transmission_dsa_results.csv", row.names = FALSE)
cat("Saved: output/csv/transmission_dsa_results.csv\n")

# Print sequential ICERs for frontier strategies across contact values
cat("\nTransmission contacts DSA — sequential ICERs for Parallel Sx+QFT (Ultra):\n")
cat(sprintf("  %-30s %12s %22s %15s\n",
    "Contacts/case/yr", "Implied beta", "Parallel Sx+QFT cost/pp", "Seq. ICER"))
for (nc in tx_dsa_contacts) {
  row <- tx_dsa_df[tx_dsa_df$n_contacts == nc &
                   tx_dsa_df$strategy == "Parallel Sx+QFT (Ultra)", ]
  cat(sprintf("  %-30s %12s %22s %15s\n",
    sprintf("%d contacts%s", nc, if (nc == n_contacts_base) " (base)" else ""),
    sprintf("\u03b2 = %.3f", nc * p_tx_per_contact),
    paste0("\u00a3", format(round(row$cost / n_c), big.mark = ",")),
    if (length(row$sequential_icer) == 0 || is.na(row$sequential_icer))
      "\u2014" else paste0("\u00a3", format(round(row$sequential_icer), big.mark = ","))))
}
cat("\n")

# Grouped bar chart: Cough+CXR + Parallel Cough+QFT + Parallel Sx+QFT × 3 contact values
tx_dsa_strats <- c("Cough+CXR (TB sx)", "Parallel Cough+QFT (Ultra)", "Parallel Sx+QFT (Ultra)")
tx_dsa_cols <- c(
  "Cough+CXR (TB sx)"          = "#c84b6a",
  "Parallel Cough+QFT (Ultra)" = "#9b5fc0",
  "Parallel Sx+QFT (Ultra)"    = "#5a1a7a"
)
tx_dsa_short <- c(
  "Cough+CXR (TB sx)"          = "Cough+CXR (TB sx)",
  "Parallel Cough+QFT (Ultra)" = "Parallel Cough+QFT (Ultra)",
  "Parallel Sx+QFT (Ultra)"    = "Parallel Sx+QFT (Ultra)"
)
tx_dsa_xlabs <- c(
  "5"  = "5 contacts/yr\n(base case)",
  "10" = "10 contacts/yr\n(mid)",
  "15" = "15 contacts/yr\n(high)"
)

tx_dsa_plot_df <- tx_dsa_df %>%
  filter(strategy %in% tx_dsa_strats) %>%
  mutate(
    strategy   = factor(strategy, levels = tx_dsa_strats),
    dominated  = dominance == "extendedly dominated",
    seq_icer_k = ifelse(!dominated, sequential_icer / 1000, NA_real_)
  )

p_tx_dsa <- ggplot(tx_dsa_plot_df,
                   aes(x = factor(n_contacts), y = seq_icer_k, fill = strategy)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, na.rm = TRUE) +
  geom_hline(yintercept = 25, linetype = "dashed", colour = "grey30", linewidth = 0.6) +
  annotate("text", x = 0.55, y = 25.4, label = "\u00a325,000/QALY (NICE WTP)",
           hjust = 0, size = 3, colour = "grey30") +
  # "Ext. dom." label — keep all 3 groups so dodge positions correctly
  geom_text(aes(y    = ifelse(dominated, 0.4, NA_real_),
                label = ifelse(dominated, "Ext. dom.", NA_character_),
                group = strategy),
            position = position_dodge(width = 0.72),
            size = 2.8, colour = "grey50", fontface = "italic",
            vjust = 0, na.rm = TRUE) +
  # Value labels on bars
  geom_text(aes(y     = seq_icer_k,
                label = ifelse(!dominated,
                               paste0("\u00a3", format(round(sequential_icer), big.mark = ",", trim = TRUE)),
                               NA_character_),
                group = strategy),
            position = position_dodge(width = 0.72),
            vjust = -0.4, size = 2.8, colour = "grey20", na.rm = TRUE) +
  scale_fill_manual(values = tx_dsa_cols, labels = tx_dsa_short, name = NULL) +
  scale_x_discrete(labels = tx_dsa_xlabs) +
  scale_y_continuous(
    labels = function(x) paste0("\u00a3", format(round(x * 1000), big.mark = ",")),
    limits = c(0, 9),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal(base_size = 12) +
  labs(
    x        = "Close contacts per untreated active TB case per year",
    y        = "Sequential ICER (\u00a3/QALY)",
    title    = NULL,
    subtitle = NULL
  ) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "bottom",
    legend.text        = element_text(size = 10)
  )

ggsave("output/transmission_dsa_bar.png", p_tx_dsa, width = 9, height = 6, dpi = 300)
cat("Saved: output/transmission_dsa_bar.png\n")

# Cost sensitivity: base (£6,055) vs upper (£12,110) — Green 2025 ERJ Open Res DSA range
cat("\nCost per secondary TB case sensitivity (5 contacts base; Parallel Sx+QFT (Ultra)):\n")
cat(sprintf("  %-18s %30s %15s\n", "Cost/case", "Parallel Sx+QFT cost/pp", "Seq. ICER"))
for (cost_test in c(cost_per_secondary_tb, as.numeric(config_list["cost_secondary_tb_upper"]))) {
  sav_test <- atb_cases_prevent * beta_tx_base * cost_test
  sr_cost <- lapply(strategy_results, function(x) {
    x$total_cost <- x$total_cost - sav_test[x$strategy_name]; x })
  tbl_cost <- calculate_icer(sr_cost)
  row_u    <- tbl_cost[tbl_cost$strategy == "Parallel Sx+QFT (Ultra)", ]
  cat(sprintf("  \u00a3%-17s %30s %15s\n",
    format(cost_test, big.mark = ","),
    paste0("\u00a3", format(round(row_u$cost / n_c), big.mark = ",")),
    if (length(row_u$sequential_icer) == 0 || is.na(row_u$sequential_icer)) "\u2014"
    else paste0("\u00a3", format(round(row_u$sequential_icer), big.mark = ","))))
}
cat("\n")

cat("\n================================================================================\n")

# =============================================================================
# PROBABILISTIC SENSITIVITY ANALYSIS (PSA)
#
# 1,000 MC simulations; params drawn from config.csv distributions; fixed seed.
# Outputs: CE plane, CEAC, ICER 95% CrIs.
# =============================================================================

set.seed(42)  # fixed seed ensures results are reproducible across runs
psa_df <- run_psa(strategies_base = strategies, n_sim = 1000)

# Save PSA raw results for supplementary materials and audit trail
write.csv(psa_df, "output/csv/psa_results.csv", row.names = FALSE)

# =============================================================================
# TRANSMISSION PSA: PROBABILISTIC SENSITIVITY ANALYSIS WITH CONTACTS
#
# n_contacts ~ Uniform(5, 15) per untreated active TB case per year.
# Implied beta_i = n_contacts_i × p_tx_per_contact (0.041), covering the
# full range from 5 (conservative) to 15 (high) contacts per case.
# This corresponds to beta ~ Uniform(0.21, 0.62), consistent with and
# slightly wider than Brooks-Pollock 2020 95% CrI (0.30-0.60).
# Active TB cases prevented held fixed from base case run; only contacts
# vary per simulation. Merged with main PSA cost draws.
# =============================================================================
cat("\nRunning transmission PSA (1,000 simulations, contacts ~ Uniform(5, 15))...\n")

# No separate seed — transmission PSA continues from RNG state after main PSA.
# Both analyses are fully reproducible from the single seed set before run_psa().
psa_tx_results <- list()
for (i in 1:1000) {
  if (i %% 200 == 0) cat(sprintf("  Transmission PSA simulation %d/1000\n", i))
  n_contacts_i <- runif(1, min = n_contacts_base, max = n_contacts_dsa_high)
  beta_i       <- n_contacts_i * p_tx_per_contact
  savings_i    <- atb_cases_prevent * beta_i * cost_per_secondary_tb
  psa_tx_results[[i]] <- tibble(
    sim        = i,
    strategy   = names(atb_cases_prevent),
    n_contacts = n_contacts_i,
    beta_tx    = beta_i,
    savings    = as.numeric(savings_i)
  )
}
psa_tx_savings <- bind_rows(psa_tx_results)

psa_tx_df <- psa_df %>%
  left_join(psa_tx_savings, by = c("sim", "strategy")) %>%
  mutate(cost_tx = cost - savings)

ref_tx <- psa_tx_df %>% filter(strategy == "No screening") %>%
  select(sim, ref_cost_tx = cost_tx, ref_qaly = qaly)

psa_tx_icer <- psa_tx_df %>%
  filter(strategy != "No screening") %>%
  left_join(ref_tx, by = "sim") %>%
  mutate(
    inc_cost_tx = cost_tx - ref_cost_tx,
    inc_qaly    = qaly - ref_qaly,
    icer_tx     = ifelse(inc_qaly > 0, inc_cost_tx / inc_qaly, NA_real_),
    ce_25k      = icer_tx < 25000,
    ce_35k      = icer_tx < 35000
  )

frontier_tx_summary <- psa_tx_icer %>%
  filter(strategy %in% c("Cough+CXR (TB sx)",
                          "Parallel Cough+QFT (Ultra)",
                          "Parallel Sx+QFT (Ultra)")) %>%
  group_by(strategy) %>%
  summarise(
    median_icer = median(icer_tx, na.rm = TRUE),
    lo_icer     = quantile(icer_tx, 0.025, na.rm = TRUE),
    hi_icer     = quantile(icer_tx, 0.975, na.rm = TRUE),
    prob_ce_25k = mean(ce_25k, na.rm = TRUE),
    prob_ce_35k = mean(ce_35k, na.rm = TRUE),
    .groups     = "drop"
  )

cat("\nTransmission PSA — frontier strategies (contacts ~ Uniform(5, 15); implied beta ~ Uniform(0.21, 0.62)):\n")
cat(sprintf("  %-22s %12s %24s %12s %12s\n",
    "Strategy", "Median ICER", "95% CrI", "P(CE \u00a325k)", "P(CE \u00a335k)"))
cat(strrep("-", 88), "\n")
for (i in seq_len(nrow(frontier_tx_summary))) {
  r <- frontier_tx_summary[i, ]
  cat(sprintf("  %-22s %12s %24s %11.0f%% %11.0f%%\n",
    r$strategy,
    paste0("\u00a3", format(round(r$median_icer), big.mark = ",")),
    paste0("\u00a3", format(round(r$lo_icer), big.mark = ","), " to \u00a3",
           format(round(r$hi_icer), big.mark = ",")),
    r$prob_ce_25k * 100,
    r$prob_ce_35k * 100))
}
write.csv(frontier_tx_summary, "output/csv/transmission_psa_summary.csv", row.names = FALSE)
cat("Saved: output/csv/transmission_psa_summary.csv\n")
cat("\n================================================================================\n")

# Summarise mean and standard deviation of cost and QALYs across simulations
psa_summary <- psa_df %>%
  group_by(strategy) %>%
  summarise(
    mean_cost = mean(cost),
    sd_cost = sd(cost),
    mean_qaly = mean(qaly),
    sd_qaly = sd(qaly),
    .groups = "drop"
  )

cat("\n")
cat("================================================================================\n")
cat("                    PSA RESULTS (1,000 simulations)                             \n")
cat("================================================================================\n\n")
cat(sprintf("%-20s %14s %14s %14s %14s\n",
            "Strategy", "Mean Cost", "SD Cost", "Mean QALYs", "SD QALYs"))
cat(sprintf("%-20s %14s %14s %14s %14s\n",
            "--------------------", "--------------", "--------------", "--------------", "--------------"))
for (i in 1:nrow(psa_summary)) {
  cat(sprintf("%-20s %13s %13s %14.1f %14.1f\n",
              psa_summary$strategy[i],
              paste0("£", format(round(psa_summary$mean_cost[i]), big.mark = ",")),
              paste0("£", format(round(psa_summary$sd_cost[i]), big.mark = ",")),
              psa_summary$mean_qaly[i],
              psa_summary$sd_qaly[i]))
}
cat("\n================================================================================\n")

# -------------------- ICER credible intervals from PSA ------------------------
# ICER 95% CrI from 2.5th/97.5th percentiles; P(CE) = proportion of sims with ICER < £25k.
ref_psa <- psa_df %>% filter(strategy == "No screening") %>%
  select(sim, ref_cost = cost, ref_qaly = qaly)

icer_psa <- psa_df %>%
  filter(strategy != "No screening") %>%
  left_join(ref_psa, by = "sim") %>%
  mutate(
    inc_cost = cost - ref_cost,
    inc_qaly = qaly - ref_qaly,
    icer = ifelse(inc_qaly != 0, inc_cost / inc_qaly, NA)
  )

icer_ci <- icer_psa %>%
  group_by(strategy) %>%
  summarise(
    mean_icer = mean(icer, na.rm = TRUE),
    median_icer = median(icer, na.rm = TRUE),
    icer_lo = quantile(icer, 0.025, na.rm = TRUE),
    icer_hi = quantile(icer, 0.975, na.rm = TRUE),
    mean_inc_cost = mean(inc_cost),
    mean_inc_qaly = mean(inc_qaly),
    prob_ce_25k = mean(inc_qaly > 0 & (inc_cost / inc_qaly) < 25000, na.rm = TRUE),
    prob_ce_35k = mean(inc_qaly > 0 & (inc_cost / inc_qaly) < 35000, na.rm = TRUE),
    prob_dominated = mean(inc_cost > 0 & inc_qaly <= 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_inc_cost)

cat("\n")
cat("================================================================================\n")
cat("                    ICER WITH 95% CREDIBLE INTERVALS (PSA)                     \n")
cat("================================================================================\n\n")
cat(sprintf("%-20s %12s %12s %24s %10s %10s %8s\n",
            "Strategy", "Mean ICER", "Median ICER", "95% CrI", "P(CE £25k)", "P(CE £35k)", "P(Dom)"))
cat(sprintf("%-20s %12s %12s %24s %10s %10s %8s\n",
            "--------------------", "------------", "------------",
            "------------------------", "----------", "----------", "--------"))
for (i in 1:nrow(icer_ci)) {
  mean_str <- paste0("£", format(round(icer_ci$mean_icer[i]), big.mark = ","))
  med_str  <- paste0("£", format(round(icer_ci$median_icer[i]), big.mark = ","))
  lo_str   <- format(round(icer_ci$icer_lo[i]), big.mark = ",")
  hi_str   <- format(round(icer_ci$icer_hi[i]), big.mark = ",")
  cri_str  <- paste0("(£", lo_str, " to £", hi_str, ")")
  cat(sprintf("%-20s %12s %12s %24s %7.1f%% %7.1f%% %7.1f%%\n",
              icer_ci$strategy[i], mean_str, med_str, cri_str,
              100 * icer_ci$prob_ce_25k[i],
              100 * icer_ci$prob_ce_35k[i],
              100 * icer_ci$prob_dominated[i]))
}
cat("\nP(CE) = Probability cost-effective at £25,000/QALY WTP threshold")
cat("\nP(Dom) = Probability dominated (more costly and fewer QALYs vs No Screening)\n")
cat("\n================================================================================\n")

# =============================================================================
# PSA SUMMARY — EFFICIENT FRONTIER STRATEGIES
#
# Extracts PSA credible intervals for the five efficient-frontier strategies
# (No Screening reference + four non-dominated strategies). Placed after the
# full 43-strategy table to provide a reader-friendly focal summary.
#
# Key contrast highlighted here:
#   Cough+CXR (TB sx):      ~100% P(CE at £25k) — ICER driven by active TB detection
#   Symptom screen+CXR:     ~100% P(CE at £25k) — similar driver
#   Parallel IGRA strategies: lower P(CE) — driven by LTBI treatment cascade
#     (initiation, completion, prevention of reactivation — all sampled in PSA)
# =============================================================================
cat("\n================================================================================\n")
cat("          PSA SUMMARY — EFFICIENT FRONTIER STRATEGIES ONLY                      \n")
cat("================================================================================\n\n")
cat("  (Full 43-strategy table above; this block highlights the narrative contrast.)\n\n")

frontier_psa <- icer_ci %>%
  filter(strategy %in% non_dominated_names) %>%
  arrange(mean_inc_cost)

cat(sprintf("  %-34s  %12s  %26s  %10s\n",
            "Strategy", "Median ICER", "95% CrI", "P(CE £25k)"))
cat(strrep("-", 88), "\n")
for (i in seq_len(nrow(frontier_psa))) {
  med_str <- if (is.na(frontier_psa$median_icer[i])) "Reference"
             else paste0("\u00a3", format(round(frontier_psa$median_icer[i]), big.mark = ","))
  lo_str  <- format(round(frontier_psa$icer_lo[i]), big.mark = ",")
  hi_str  <- format(round(frontier_psa$icer_hi[i]), big.mark = ",")
  cri_str <- paste0("(\u00a3", lo_str, " to \u00a3", hi_str, ")")
  pce_str <- if (is.na(frontier_psa$prob_ce_25k[i])) "—"
             else sprintf("%.1f%%", 100 * frontier_psa$prob_ce_25k[i])
  cat(sprintf("  %-34s  %12s  %26s  %10s\n",
              frontier_psa$strategy[i], med_str, cri_str, pce_str))
}

# Derive summary stats dynamically from data for interpretation block
.igra_row   <- frontier_psa[frontier_psa$strategy == "Parallel Sx+QFT (Ultra)", ]
.igra_pce   <- if (nrow(.igra_row) > 0) sprintf("%.1f%%", 100 * .igra_row$prob_ce_25k) else "—"
.igra_pdom  <- if (nrow(.igra_row) > 0) sprintf("%.1f%%", 100 * .igra_row$prob_dominated) else "—"

cat("\n")
cat("  Interpretation:\n")
cat("  - Non-IGRA frontier strategies (Cough+CXR, Symptom screen+CXR) have near-certain\n")
cat("    CE — ICER driven by active TB detection, which is stable across PSA draws.\n")
cat(sprintf("  - Parallel Sx+QFT (Ultra) P(CE £25k) = %s; P(dominated vs No Screening) = %s.\n",
            .igra_pce, .igra_pdom))
cat("    Wider uncertainty because QALY advantage depends on LTBI treatment cascade\n")
cat("    (initiation, completion, prevention of reactivation), each sampled independently.\n")
cat("\n================================================================================\n")

# Save ICER credible interval table for results reporting
icer_ci_out <- icer_ci %>%
  mutate(
    mean_inc_cost_pp = mean_inc_cost / n_c,
    mean_inc_qaly_pp = mean_inc_qaly / n_c
  )
write.csv(icer_ci_out, "output/csv/icer_confidence_intervals.csv", row.names = FALSE)
cat("ICER CI table saved to output/csv/icer_confidence_intervals.csv\n")

# -------------------- ICER forest plot with 95% credible intervals ------------
# ICERs are capped at ±£200,000 for axis readability; strategies with very
# high uncertainty (e.g. dominated in most simulations) can produce extreme
# ICER values that would otherwise compress the scale. Capping is standard
# practice in health economic visualisation and does not affect interpretation.
icer_ci_plot <- icer_ci %>%
  mutate(
    median_icer_cap = pmax(pmin(median_icer, 200000), -200000),
    icer_lo_cap     = pmax(pmin(icer_lo, 200000), -200000),
    icer_hi_cap     = pmax(pmin(icer_hi, 200000), -200000),
    strategy        = reorder(strategy, median_icer)
  )

p_icer_forest <- ggplot(icer_ci_plot, aes(x = median_icer_cap, y = strategy)) +
  geom_vline(xintercept = 25000, linetype = "dashed", color = "#7a1a4a", linewidth = 0.8) +
  geom_vline(xintercept = 35000, linetype = "dashed", color = "#c84b6a", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(xmin = icer_lo_cap, xmax = icer_hi_cap),
                width = 0.3, linewidth = 0.8, color = "grey40", orientation = "y") +
  geom_point(size = 4, color = "#c84b6a") +
  annotate("label", x = 25000, y = Inf, vjust = 1.3, label = "NICE £25k\nWTP threshold",
           color = "#7a1a4a", fill = "white", alpha = 0.85,
           size = 3, hjust = 1.05, fontface = "bold", linewidth = 0.3) +
  annotate("label", x = 35000, y = Inf, vjust = 1.3, label = "NICE £35k\nWTP threshold",
           color = "#c84b6a", fill = "white", alpha = 0.85,
           size = 3, hjust = -0.05, fontface = "bold", linewidth = 0.3) +
  scale_x_continuous(labels = function(x) paste0("£", trimws(format(x, big.mark = ",", scientific = FALSE)))) +
  theme_minimal(base_size = 14) +
  labs(x = "ICER (£/QALY)", y = "",
       title = NULL,
       subtitle = "Median ICER from 1,000 PSA simulations vs No Screening") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9))

ggsave("output/icer_forest_plot.png", p_icer_forest, width = 14, height = 14, dpi = 300)

# -------------------- Probability of cost-effectiveness stacked bar chart -----
# Each bar shows the proportion of PSA simulations in which the strategy was
# cost-effective (ICER < £25,000/QALY), not cost-effective, or dominated.
pce_data <- icer_ci %>%
  select(strategy, prob_ce_25k, prob_dominated) %>%
  mutate(
    prob_not_ce = 1 - prob_ce_25k - prob_dominated,
    strategy = reorder(strategy, prob_ce_25k)
  ) %>%
  pivot_longer(cols = c(prob_ce_25k, prob_dominated, prob_not_ce),
               names_to = "outcome", values_to = "prob") %>%
  mutate(outcome = factor(outcome,
    levels = c("prob_ce_25k", "prob_not_ce", "prob_dominated"),
    labels = c("Cost-effective (ICER < £25k)", "Not cost-effective", "Dominated")))

p_pce_bar <- ggplot(pce_data, aes(x = strategy, y = prob, fill = outcome)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Cost-effective (ICER < £25k)" = "#c84b6a",
                                "Not cost-effective" = "#e8a0b8",
                                "Dominated" = "#2d0a2e")) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "Proportion of PSA simulations",
       title = NULL,
       subtitle = "At £25,000/QALY WTP threshold | 1,000 PSA simulations",
       fill = "") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(size = 9))

ggsave("output/prob_cost_effective.png", p_pce_bar, width = 13, height = 14, dpi = 300)

# -------------------- Cost-effectiveness plane with 95% confidence ellipses ---
# Each cloud of points represents one strategy across 1,000 PSA simulations.
# The 95% ellipse captures joint cost-QALY uncertainty. Diamond markers show
# the mean position. The dotted line is the £25,000/QALY WTP frontier.
icer_psa_pp <- icer_psa %>%
  mutate(inc_cost_pp = inc_cost / n_c, inc_qaly_pp = inc_qaly / n_c)

icer_means <- icer_psa_pp %>%
  group_by(strategy) %>%
  summarise(inc_cost_pp = mean(inc_cost_pp), inc_qaly_pp = mean(inc_qaly_pp), .groups = "drop")

all_strat_names  <- unique(psa_df$strategy)
ce_plane_colours <- build_strategy_colours(all_strat_names)

# For ellipse / CE plane / CEAC: filter to non-dominated strategies only.
# 35 strategies produce unreadable overlapping clouds and legends.
# Standard HTA practice: show efficient frontier + reference on these plots;
# all 43 are shown in the strategy space scatter and forest/bar plots.
nd_strats_psa    <- non_dominated_names[non_dominated_names != "No screening"]
icer_psa_nd      <- icer_psa_pp %>% filter(strategy %in% nd_strats_psa)
icer_means_nd    <- icer_means   %>% filter(strategy %in% nd_strats_psa)
ce_ellipse_cols  <- build_strategy_colours(nd_strats_psa)

# Distinct colours for CE plane — match CEAC palette for consistency
ce_ellipse_named_cols <- c(
  "Cough+CXR (TB sx)"          = "#c84b6a",
  "Symptom screen+CXR"         = "#9b5fc0",
  "Parallel Sx+QFT (Ultra)"    = "#5a1a7a"
)
extra_ellipse <- setdiff(nd_strats_psa, names(ce_ellipse_named_cols))
if (length(extra_ellipse) > 0) {
  extra_ce_cols <- setNames(project_pal(length(extra_ellipse)), extra_ellipse)
  ce_ellipse_named_cols <- c(ce_ellipse_named_cols, extra_ce_cols)
}

p_ce_ellipse <- ggplot(icer_psa_nd, aes(x = inc_qaly_pp, y = inc_cost_pp, color = strategy)) +
  geom_point(alpha = 0.15, size = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  geom_point(data = icer_means_nd, size = 4, shape = 18) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(intercept = 0, slope = 25000, linetype = "dotted", color = "black", linewidth = 1) +
  # Annotate the Cough+CXR cluster which sits at the origin and is too small
  # to see at this scale (inc. QALYs ~0.006 pp; inc. cost ~£7 pp).
  annotate("segment", x = 0.005, xend = 0.001, y = 30, yend = 5,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
           colour = "#c84b6a", linewidth = 0.8) +
  annotate("text", x = 0.006, y = 38,
           label = "Cough+CXR (TB sx)\n(inc. cost ~\u00a37, inc. QALY ~0.006;\nsee zoomed plot)",
           colour = "#c84b6a", size = 3, hjust = 0) +
  scale_colour_manual(values = ce_ellipse_named_cols) +
  scale_x_continuous(labels = function(x) format(round(x, 4), nsmall = 4, scientific = FALSE)) +
  scale_y_continuous(labels = function(x) paste0(ifelse(x < 0, "-", ""), "\u00a3", scales::comma(abs(x)))) +
  theme_minimal(base_size = 14) +
  labs(x = "Incremental QALYs per person",
       y = "Incremental cost per person (\u00a3)",
       title = NULL,
       subtitle = "Efficient frontier strategies only | vs No Screening | Dotted = \u00a325,000/QALY | 1,000 PSA simulations",
       color = "Strategy") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(),
        legend.position = "right")

ggsave("output/ce_plane_ellipses.png", p_ce_ellipse, width = 13, height = 8, dpi = 300)

# Zoomed CE plane — shows the Cough+CXR (TB sx) cluster near the origin
# which is invisible in the full-scale ellipses plot above.
p_ce_ellipse_zoom <- ggplot(icer_psa_nd, aes(x = inc_qaly_pp, y = inc_cost_pp, color = strategy)) +
  geom_point(alpha = 0.25, size = 1) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  geom_point(data = icer_means_nd, size = 4, shape = 18) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(intercept = 0, slope = 25000, linetype = "dotted", color = "black", linewidth = 1) +
  coord_cartesian(xlim = c(-0.005, 0.020), ylim = c(-20, 60)) +
  scale_colour_manual(values = ce_ellipse_named_cols) +
  scale_x_continuous(labels = function(x) format(round(x, 4), nsmall = 4, scientific = FALSE)) +
  scale_y_continuous(labels = function(x) paste0(ifelse(x < 0, "-", ""), "\u00a3", scales::comma(abs(x)))) +
  theme_minimal(base_size = 14) +
  labs(x = "Incremental QALYs per person",
       y = "Incremental cost per person (\u00a3)",
       title = NULL,
       subtitle = "Zoomed to origin region | Parallel IGRA ellipses extend far right (see full-scale plot)",
       color = "Strategy") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(),
        legend.position = "right")

ggsave("output/ce_plane_ellipses_zoom.png", p_ce_ellipse_zoom, width = 13, height = 8, dpi = 300)
cat("Saved: output/ce_plane_ellipses_zoom.png\n")

cat("ICER visualisations saved to output/ folder\n")

# CE plane and CEAC: filter PSA data to non-dominated strategies for readability.
# "No screening" (reference) must be included for incremental calculations.
psa_df_nd <- psa_df %>% filter(strategy %in% non_dominated_names)

# Cost-effectiveness plane: scatter of incremental cost vs incremental QALY
# across all 1,000 simulations for efficient frontier strategies vs no screening.
p_ce_plane <- plot_ce_plane(psa_df_nd, reference = "No screening", wtp = 25000)

# Cost-effectiveness acceptability curve (CEAC): probability each strategy is
# optimal across WTP thresholds. Efficient frontier strategies only.
p_ceac <- plot_ceac(psa_df_nd, reference = "No screening")

ggsave("output/ce_plane.png", p_ce_plane, width = 12, height = 8, dpi = 300)
ggsave("output/ceac.png",     p_ceac,     width = 12, height = 7, dpi = 300)

cat("\nPSA plots saved to output/ folder\n")

# =============================================================================
# SUPPLEMENTARY EXCEL FILES — S4 and S5
# =============================================================================
# S4: All 43 strategies × all 3 conditions (base case + 2 scenarios)
# S5: All sensitivity analyses (DSA, PSA, structural SAs) × all scenarios
#
# Both files are written to docs/ with a date-stamped filename.
# All data come from CSVs already written to output/csv/ in this run.
# =============================================================================

library(openxlsx)

.today_str <- format(Sys.Date(), "%Y-%m-%d")

# Helper: read a CSV and write it as a formatted table into a workbook sheet.
# Returns the data frame (invisibly) so callers can inspect it if needed.
.write_supp_sheet <- function(wb, sheet_name, csv_path, description) {
  df <- tryCatch(
    read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) {
    cat(sprintf("  [SKIP] %s not found\n", csv_path))
    return(invisible(NULL))
  }
  addWorksheet(wb, sheetName = sheet_name)

  # Row 1: description note
  writeData(wb, sheet_name, description,
            startRow = 1, startCol = 1, colNames = FALSE)
  addStyle(wb, sheet_name,
           createStyle(fontSize = 9, fontColour = "#555555", wrapText = TRUE),
           rows = 1, cols = 1, gridExpand = FALSE)
  setRowHeights(wb, sheet_name, rows = 1, heights = 28)

  # Row 2+: data table with auto-filter
  writeDataTable(wb, sheet_name, df,
                 startRow = 2, startCol = 1,
                 tableStyle = "TableStyleLight9",
                 withFilter  = TRUE)
  setColWidths(wb, sheet_name,
               cols   = seq_len(ncol(df)),
               widths = "auto")
  invisible(df)
}

# ---------------------------------------------------------------------------
# S4 — All strategies, all conditions
# ---------------------------------------------------------------------------
.s4_path <- file.path("docs", "S4. Full cost-effectiveness results.xlsx")
.wb4 <- createWorkbook()
modifyBaseFont(.wb4, fontSize = 10, fontName = "Calibri")

.write_supp_sheet(
  .wb4, "Base case",
  "output/csv/icer_table_basecase.csv",
  paste0("S4a · Base case — all 43 strategies · 55-year lifetime horizon · ",
         "n = 100,000 cohort · 2023/24 GBP · ", .today_str,
         " · Columns: cost and QALYs per cohort; cost/QALY per person; ",
         "ICER vs No Screening; sequential ICER; dominance classification; ",
         "LTBI detected; active TB detected.")
)

.write_supp_sheet(
  .wb4, "Scenario 1 — Inaction",
  "output/csv/icer_table_inaction.csv",
  paste0("S4b · Scenario 1 (Cost of unmanaged TB) — all 43 strategies · ",
         "Unmanaged active TB: £150/month undiagnosed; £275/month untreated; ",
         "£175/month lost-to-follow-up · All other parameters at base-case values · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb4, "Scenario 2 — Transmission",
  "output/csv/icer_table_transmission.csv",
  paste0("S4c · Scenario 2 (Secondary cases prevention) — all 43 strategies · ",
         "5 contacts/case/year (conservative base) × 4.1% annual transmission risk = ",
         "β = 0.205 infections/year · £6,055/secondary case prevented · ",
         "All other parameters at base-case values · 2023/24 GBP · ", .today_str)
)

saveWorkbook(.wb4, .s4_path, overwrite = TRUE)
cat(sprintf("Saved: %s\n", .s4_path))

# ---------------------------------------------------------------------------
# S5 — All sensitivity analyses, DSA and PSA
# ---------------------------------------------------------------------------
.s5_path <- file.path("docs", "S5. Sensitivity analyses.xlsx")
.wb5 <- createWorkbook()
modifyBaseFont(.wb5, fontSize = 10, fontName = "Calibri")

# --- PSA ---
.write_supp_sheet(
  .wb5, "PSA — Base case summary",
  "output/csv/icer_confidence_intervals.csv",
  paste0("S5a · PSA base case — ICER 95% credible intervals and P(CE) for all 43 strategies · ",
         "1,000 simulations · Beta distributions for proportions/utilities; ",
         "Gamma distributions for costs · WTP thresholds £25,000 and £35,000/QALY · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "PSA — Raw simulations",
  "output/csv/psa_results.csv",
  paste0("S5b · PSA raw results — cost and QALYs per cohort (100,000) for each of ",
         "1,000 simulations × 43 strategies (43,000 rows) · ",
         "sim = simulation index; cost = total cohort cost (£); qaly = total cohort QALYs · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "PSA — Transmission scenario",
  "output/csv/transmission_psa_summary.csv",
  paste0("S5c · Scenario 2 (Transmission) PSA summary — frontier strategies · ",
         "Contacts per case/year drawn from Uniform(5, 15) across 1,000 simulations · ",
         "β = contacts × 4.1% annual transmission risk · £6,055/secondary case · ",
         "2023/24 GBP · ", .today_str)
)

# --- DSA ---
.write_supp_sheet(
  .wb5, "DSA — Cough+CXR (TB sx)",
  "output/csv/dsa_results_coughcxr.csv",
  paste0("S5d · One-way DSA — Cough+CXR (TB sx) vs No Screening · ",
         "NMB at £25,000/QALY WTP · Each parameter varied independently between ",
         "its lower and upper bound; all others held at base-case values · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "DSA — Symptom screen+CXR",
  "output/csv/dsa_results_symscrCXR.csv",
  paste0("S5e · One-way DSA — Symptom screen+CXR vs No Screening · ",
         "NMB at £25,000/QALY WTP · Each parameter varied independently between ",
         "its lower and upper bound; all others held at base-case values · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "DSA — Parallel Sx+QFT (Ultra)",
  "output/csv/dsa_results.csv",
  paste0("S5f · One-way DSA — Parallel Sx+QFT (Ultra) vs No Screening · ",
         "NMB at £25,000/QALY WTP · Each parameter varied independently between ",
         "its lower and upper bound; all others held at base-case values · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "DSA — Transmission contacts",
  "output/csv/transmission_dsa_results.csv",
  paste0("S5g · Scenario 2 DSA — close contacts per active TB case · ",
         "Values: 3, 5 (base), 7, 10, 15 contacts/case/year · ",
         "β = contacts × 4.1% per year · £6,055/secondary case · ",
         "2023/24 GBP · ", .today_str)
)

# --- Structural SAs ---
.write_supp_sheet(
  .wb5, "SA — LTBI completion",
  "output/csv/ltbi_completion_sensitivity.csv",
  paste0("S5h · Structural SA — LTBI treatment completion · ",
         "Scenarios: trial 76.7% (base; UKHSA 2024 TB Prevention report) vs ",
         "real-world 55.5% (UKHSA 2024) · 2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "SA — LTBI prevalence",
  "output/csv/ltbi_prevalence_sensitivity.csv",
  paste0("S5i · Structural SA — LTBI prevalence at entry · ",
         "Scenarios: base 17.8% (UKHSA 2021 migrant cohort) vs ",
         "15.1% (UKHSA 2025 updated estimate) · 2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "SA — Active TB prevalence",
  "output/csv/active_tb_prev_sensitivity.csv",
  paste0("S5j · Structural SA — Active TB prevalence at entry · ",
         "Scenarios: base 1.0%, high-risk 0.44%, pooled low-burden 0.215% · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "SA — IGRA uptake",
  "output/csv/igra_uptake_sensitivity.csv",
  paste0("S5k · Structural SA — IGRA programme uptake rate · ",
         "Varied 70% to 100% · Applies to all parallel IGRA strategies · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "SA — IGRA specificity",
  "output/csv/igra_specificity_sensitivity.csv",
  paste0("S5l · Structural SA — IGRA test specificity · ",
         "QFT-GIT base 0.96; T-SPOT.TB base 0.93; varied 0.90 to 1.00 · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "SA — IGRA programme cost",
  "output/csv/igra_programme_cost_sensitivity.csv",
  paste0("S5m · Structural SA — IGRA programme delivery cost per person screened · ",
         "Varied across plausible range (£4 to £60) · ",
         "Cost added to initial screening cost only; does not affect Markov transitions · ",
         "2023/24 GBP · ", .today_str)
)

saveWorkbook(.wb5, .s5_path, overwrite = TRUE)
cat(sprintf("Saved: %s\n", .s5_path))

cat("\n================================================================================\n")
cat("                    ALL ANALYSES COMPLETE                                       \n")
cat("================================================================================\n")
