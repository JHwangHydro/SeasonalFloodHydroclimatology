#===========================================================
# 1. Fit one subset model and calculate station-level R2
#===========================================================
get_R2_by_station <- function(fixed_predictors, data) {
  
  if (length(fixed_predictors) == 0) {
    
    # empty model (null): predict mean by group
    
    formula_obj <- qmax ~ 1 + (1 | group_huc2) + (1 | group_huc4) + (1 | group_station)
    
  } else {
    
    fixed_predictors <- unique(fixed_predictors)
    
    fixed_str <- paste(fixed_predictors, collapse = " + ")
    
    if ("qpre" %in% fixed_predictors){
      station_term <- "(qpre || group_station)"
    }else{
      station_term <- "(1 | group_station)"
    }
    
    formula_obj <- as.formula(
      paste0(
        "qmax ~ ",
        fixed_str,
        " + (",
        fixed_str,
        " || group_huc2) + (",
        fixed_str,
        " || group_huc4) + ",
        station_term
      )
    )
    
  }
  
  # Run model
  mod <- lmer(
    formula_obj,
    data = data,
    REML = F,
    control = lmerControl(
      optCtrl = list(maxfun = 2e5),
      optimizer = "bobyqa",
    )
  )
  
  # Calculate R2
  predicted <- predict(mod, newdata = data)
  
  result_df <- data %>%
    transmute(
      group_station = group_station,
      qmax = qmax,
      pre = predicted
    )
  
  r2_by_station <- result_df %>%
    dplyr::group_by(group_station) %>%
    dplyr::summarise(
      SSE = sum((qmax - pre)^2),
      SST = sum((qmax - mean(qmax))^2),
      R2 = if_else(
        is.finite(SST) & SST > 0,
        1 - SSE / SST,
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::select(group_station, R2)
  
  # Enforce fixed station order
  stations <- sort(unique(data$group_station))
  
  r2_by_station <- tibble(
    group_station = stations
  ) %>%
    left_join(
      r2_by_station,
      by = "group_station"
    )
  
  # Model diagnostic
  singular_fit <- isSingular(mod, tol = 1e-4)
  
  convergence_message <- mod@optinfo$conv$lme4$messages
  if (!is.null(convergence_message)) {
    # Exclude singular-fit messages
    convergence_message <- convergence_message[
      !grepl("singular", convergence_message, ignore.case = TRUE)
    ]
    if (length(convergence_message) == 0) {
      convergence_message <- NA_character_
    } else {
      convergence_message <- paste(
        convergence_message,
        collapse = "; "
      )
    }
  } else {
    convergence_message <- NA_character_
  }
  
  converged <- is.na(convergence_message)
  
  # Output
  list(
    R2 = r2_by_station,
    model = mod,
    success = converged && all(is.finite(r2_by_station$R2)),
    singular = singular_fit,
    convergence_message = convergence_message
  )
  
}


#===========================================================
# 2. General dominance analysis
#===========================================================
DA_function <- function(
    data,
    groups_list
) {
  
  # Input checks
  group_names <- names(groups_list)
  n_groups <- length(group_names)
  
  all_predictors <- unique(unlist(groups_list))
  
  required_vars <- unique(
    c(
      "qmax",
      "group_huc2",
      "group_huc4",
      "group_station",
      all_predictors
    )
  )
  
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(
      "The following variables are missing from data: ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  
  # Use exactly the same observations for every model
  data_da <- data %>%
    dplyr::select(all_of(required_vars)) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::mutate(
      group_huc2 = droplevels(factor(group_huc2)),
      group_huc4 = droplevels(factor(group_huc4)),
      group_station = droplevels(factor(group_station))
    )
  
  stations <- sort(unique(data_da$group_station))
  n_station <- length(stations)
  
  
  # Cache models so identical subsets are not fitted repeatedly
  model_cache <- new.env(parent = emptyenv())
  
  predictor_key <- function(predictors) {
    
    predictors <- sort(unique(predictors))
    
    if (length(predictors) == 0) {
      return("EMPTY")
    }
    
    paste(predictors, collapse = " + ")
  }
  
  get_cached_R2 <- function(predictors) {
    
    key <- predictor_key(predictors)
    
    if (exists(key, envir = model_cache, inherits = FALSE)) {
      return(get(key, envir = model_cache))
    }
    
    fit_result <- get_R2_by_station(
      fixed_predictors = predictors,
      data = data_da
    )
    
    assign(
      key,
      fit_result,
      envir = model_cache
    )
    
    fit_result
  }
  
  
  # Fit null model
  empty_fit <- get_cached_R2(character(0))
  
  if (!empty_fit$success) {
    stop("The empty model failed to fit or produce valid station-level R2.")
  }
  
  empty_r2 <- empty_fit$R2
  
  
  # Storage
  addR2_results_by_size <- setNames(
    vector("list", n_groups),
    group_names
  )
  
  for (g in group_names) {
    addR2_results_by_size[[g]] <- list()
  }
  
  diagnostics <- list()
  
  raw_deltaR2 <- setNames(
    vector("list", n_groups),
    group_names
  )
  
  for (g in group_names) {
    raw_deltaR2[[g]] <- list()
  }
  
  
  # Helper for one nested-model comparison
  compare_models <- function(
    predictors_with,
    predictors_without,
    comparison_name
  ) {
    
    fit_with <- get_cached_R2(predictors_with)
    fit_without <- get_cached_R2(predictors_without)
    
    valid_with <- fit_with$success
    valid_without <- fit_without$success
    
    diagnostics[[length(diagnostics) + 1]] <<- tibble(
      comparison = comparison_name,
      with_predictors = predictor_key(predictors_with),
      without_predictors = predictor_key(predictors_without),
      with_success = fit_with$success,
      without_success = fit_without$success,
      with_singular = fit_with$singular,
      without_singular = fit_without$singular,
      with_convergence_message = fit_with$convergence_message,
      without_convergence_message = fit_without$convergence_message
    )
    
    if (!valid_with || !valid_without) {
      return(
        tibble(
          group_station = stations,
          deltaR2 = NA_real_
        )
      )
    }
    
    tibble(
      group_station = stations
    ) %>%
      left_join(
        fit_with$R2 %>%
          rename(R2_with = R2),
        by = "group_station"
      ) %>%
      left_join(
        fit_without$R2 %>%
          rename(R2_without = R2),
        by = "group_station"
      ) %>%
      mutate(
        # Do not truncate negative increments
        deltaR2 = R2_with - R2_without
      ) %>%
      dplyr::select(group_station, deltaR2)
  }
  
  
  # k = 0: group alone versus empty model
  for (g in group_names) {
    
    comparison <- compare_models(
      predictors_with = groups_list[[g]],
      predictors_without = character(0),
      comparison_name = paste0(g, " | empty")
    )
    
    raw_deltaR2[[g]][["k0"]] <- comparison$deltaR2
    
    deltaR2_k0 <- comparison$deltaR2
    deltaR2_k0[deltaR2_k0 < 0] <- 0
    
    addR2_results_by_size[[g]][["k0"]] <- deltaR2_k0
    
  }
  
  
  # k = 1, ..., n_groups - 1
  if (n_groups > 1) {
    
    for (k in seq_len(n_groups - 1)) {
      
      for (g in group_names) {
        
        other_groups <- setdiff(group_names, g)
        
        combinations <- combn(
          other_groups,
          k,
          simplify = FALSE
        )
        
        deltaR2_mat <- matrix(
          NA_real_,
          nrow = n_station,
          ncol = length(combinations),
          dimnames = list(
            as.character(stations),
            vapply(
              combinations,
              function(x) paste(x, collapse = " + "),
              character(1)
            )
          )
        )
        
        for (j in seq_along(combinations)) {
          
          combo <- combinations[[j]]
          
          predictors_without <- unique(
            unlist(
              groups_list[combo],
              use.names = FALSE
            )
          )
          
          predictors_with <- unique(
            c(
              predictors_without,
              groups_list[[g]]
            )
          )
          
          comparison_name <- paste0(
            g,
            " | {",
            paste(combo, collapse = ", "),
            "}"
          )
          
          comparison <- compare_models(
            predictors_with = predictors_with,
            predictors_without = predictors_without,
            comparison_name = comparison_name
          )
          
          deltaR2_mat[, j] <- comparison$deltaR2
        }
        
        raw_deltaR2[[g]][[paste0("k", k)]] <- deltaR2_mat
        
        deltaR2_mat[deltaR2_mat < 0] <- 0
        
        # A valid general-dominance result requires every
        # subset comparison for this subset size.
        if (anyNA(deltaR2_mat)) {
          
          deltaR2_k <- rep(
            NA_real_,
            n_station
          )
          
        } else {
          
          deltaR2_k <- rowMeans(deltaR2_mat)
        }
        
        addR2_results_by_size[[g]][[
          paste0("k", k)
        ]] <- deltaR2_k
      }
    }
  }
  
  
  # Average equally over subset sizes
  dominance_by_station <- tibble(
    group_station = stations
  )
  
  for (g in group_names) {
    
    component_matrix <- do.call(
      rbind,
      addR2_results_by_size[[g]]
    )
    
    dominance_value <- if (anyNA(component_matrix)){
      rep(NA_real_, n_station)
    }else{
      colMeans(component_matrix)
    }
    
    dominance_by_station[[g]] <- dominance_value
  }
  
  
  dominance_percentage <- dominance_by_station %>%
    mutate(
      total_dominance = rowSums(
        across(-group_station),
        na.rm = TRUE
      )
    ) %>%
    mutate(
      across(
        -c(group_station, total_dominance),
        ~ .x / total_dominance * 100
      )
    ) %>%
    dplyr::select(-total_dominance)
  
  
  # Model and comparison diagnostics
  diagnostics_df <- if (length(diagnostics) > 0) {
    bind_rows(diagnostics)
  } else {
    tibble()
  }
  
  cached_model_names <- ls(
    envir = model_cache,
    all.names = TRUE
  )
  
  model_diagnostics <- map_dfr(
    cached_model_names,
    function(key) {
      
      fit <- get(
        key,
        envir = model_cache
      )
      
      tibble(
        model = key,
        success = fit$success,
        singular = fit$singular,
        convergence_message = fit$convergence_message
      )
    }
  )
  
  
  # Return
  list(
    dominance_by_station = dominance_by_station,
    dominance_percentage = dominance_percentage,
    increments_by_subset_size = addR2_results_by_size,
    raw_increments = raw_deltaR2,
    model_diagnostics = model_diagnostics,
    comparison_diagnostics = diagnostics_df,
    analysis_data = data_da
  )
}
