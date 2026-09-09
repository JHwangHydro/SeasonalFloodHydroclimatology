setwd()
df_winter <- read.csv("df_winter.csv", header = T)
# qmax: monthly 3-day maxima
# qpre: monthly mean flow from precedent month
# ENSO: Nino 3.4 (3-month rolling average)
# AMO: AMO (3-month rolling average)
# NAO: NAO (3-month rolling average)
# PDO: PDO (3-month rolling average)
# PNA: PNA (3-month rolling average)
# MJO1: First principal component of rMII (3-month rolling average)
# MJO2: Second principal component of rMII 3.4 (3-month rolling average)
# df is nested for HUC2 (group_huc2), HUC4 (group_huc4), and station (group_station)

# Base model ----
mod0 <- lmer(
  qmax ~ qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 +
    (qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 || group_huc2) +
    (1 | group_huc4) +
    (1 | group_station),
  data = df,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

# Station model ----
mod_station <- lmer(
  qmax ~ qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 +
    (qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 || group_huc2) +
    (1 | group_huc4) +
    (qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 || group_station),
  data = df,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

# HUC4 model ----
mod_huc4 <- lmer(
  qmax ~ qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 +
    (qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 || group_huc2) +
    (qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 || group_huc4) +
    (1 | group_station),
  data = df,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

# Full model ----
mod_full <- lmer(
  qmax ~ qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 +
    (qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 || group_huc2) +
    (qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 || group_huc4) +
    (qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 || group_station),
  data = df,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

# Model Comparison ----
data.frame(
  model = c("mod0", "mod_station", "mod_huc4", "mod_full"),
  AIC = c(AIC(mod0), AIC(mod_station), AIC(mod_huc4), AIC(mod_full)),
  BIC = c(BIC(mod0), BIC(mod_station), BIC(mod_huc4), BIC(mod_full)),
  logLik = c(
    as.numeric(logLik(mod0)),
    as.numeric(logLik(mod_station)),
    as.numeric(logLik(mod_huc4)),
    as.numeric(logLik(mod_full))
  ),
  singular = c(
    isSingular(mod0),
    isSingular(mod_station),
    isSingular(mod_huc4),
    isSingular(mod_full)
  )
)


## Final model ----
mod_fin <- lmer(
  qmax ~ qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 +
    (qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 || group_huc2) +
    (qpre + ENSO + AMO + NAO + PDO + PNA + pc1 + pc2 || group_huc4) +
    (qpre || group_station),
  data = df_winter,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

## Dominance Analysis
source("DA_function.R")

group_list <- list(
  P1 = "qpre",
  P2 = c("ENSO", "PDO", "pc1", "pc2"),
  P3 = "AMO",
  P4 = c("NAO", "PNA")
)

DA_result <- DA_function(
  data = df_winter,
  groups_list = group_list
)