#' Identify milk loss events and resilience indicators from daily milk yields
#'
#' @param data A data frame containing the observed and predicted daily milking records.
#' @param id_col The name of the column containing the individual IDs.
#' @param dim_col The name of the column containing the days in milk.
#' @param MY_col The name of the column containing the observed milk yield.
#' @param MY_pred The name of the column containing the predicted milk yield (baseline).
#' @param dim_start The first day in milk to consider when identifying milk loss events and resilience indicators.
#' @param dim_end The last day in milk to consider when identifying milk loss events and resilience indicators.
#' @param rec_mode How "recovery" is defined. One of:
#'   "pctbase": recovery when the observed value reaches a given fraction of the baseline (\code{rec}), for a given number of consecutive days (\code{stick});
#'   "band": recovery when the observation is inside a tolerance band around the baseline (+/- \code{tol}), for at least \code{stick} consecutive days;
#'   "resid": recovery when the residual has improved enough from the nadir (by a fraction \code{rec} of the nadir's absolute residual) for \code{stick} consecutive days.
#' @param drop_pct Minimum relative drop from the anchor (baseline reference) to accept an episode.
#' @param min_len Minimum number of consecutive days with negative residuals required to define an episode.
#' @param tol Used when the "band" mode is selected. Half-width of the tolerance band around baseline in relative terms.
#' @param stick Minimum number of consecutive days in recovery to consider an episode finished.
#' @param rec Minimum relative recovery from the nadir to finish an episode (used in "pctbase" and "resid" modes).
#'
#' @return A list with two data frames:
#' \itemize{
#'   \item \code{episodes}: individual milk loss events and their resilience indicators;
#'   \item \code{aggregates}: milk loss events aggregated per individual.
#' }
#' The resilience indicators identified are described in the Details section.
#'
#' @details
#' The function computes several descriptors of milk-yield perturbation episodes.
#'
#' \strong{1) Nadir (day of minimum)}
#'
#' The worst day inside the episode (deepest point of the perturbation).
#'
#' \code{t_hat = argmin_{t in [t_start, t_end]} obs(t)}
#'
#' \code{Nadir = obs(t_hat)}
#'
#' where \code{t_start} and \code{t_end} are the episode boundaries.
#'
#' \strong{2) Amplitude (drop)}
#'
#' Depth of the dip relative to the baseline at the episode start.
#'
#' \code{A = baseline(t_start) - obs(t_hat)}
#'
#' Some variants use \code{baseline(t_hat)} instead of \code{baseline(t_start)}; here the
#' start of the episode is used as the reference.
#'
#' \strong{3) ML_per_event (AUD)}
#'
#' Total milk lost (in baseline units) over the episode, i.e., the integrated milk deficit.
#'
#' \code{ML_per_event = AUD = sum_{t=t_start..t_end} [baseline(t) - obs(t)]}
#'
#' In discrete data, AUD is computed with day-weighting: each observation contributes
#'
#' \code{(baseline(t) - obs(t)) * delta_days}
#'
#' where \code{delta_days} is the gap to the next observed DIM (last day weight = 1).
#'
#' \strong{4) Time-to-baseline (TTB)}
#'
#' Time after the nadir until the profile returns to (and stays near) the baseline.
#'
#' Recovery is declared when \code{obs(t)} re-enters a tolerance band around the
#' baseline and stays there for \code{stick} consecutive days (controlled by \code{tol} and \code{stick}).
#'
#' We find the smallest \code{tau >= 0} such that for all \code{u} in the interval from
#' \code{t_hat + tau} to \code{t_hat + tau + stick - 1}:
#'
#' \code{abs(obs(u) - baseline(u)) <= tol * baseline(u)}
#'
#' Then: \code{TTB = tau}.
#'
#' If this condition is never satisfied before DIM 305, \code{TTB} is set to \code{NA}
#' (right-censored).
#'
#' \strong{5) Recovery half-life (t_1_2)}
#'
#' Earliest time after nadir when half of the drop has been recovered.
#'
#' With amplitude \code{A} as above, define the half-recovery level:
#'
#' \code{L_half = baseline(t_start) - A / 2}
#'
#' Then:
#'
#' \code{t_1_2 = min{tau >= 0 : obs(t_hat + tau) >= L_half}}
#'
#' \strong{6) Slopes (decline and recovery)}
#'
#' Average daily change during the decline into the nadir and during early
#' recovery, summarizing the episode shape.
#'
#' For a \code{K}-day local window:
#'
#' \code{DeclineSlope = (obs(min(t_hat, t_start + K)) - obs(t_start)) / (min(t_hat, t_start + K) - t_start)}
#'
#' \code{RecoverySlope = (obs(min(t_end, t_hat + K)) - obs(t_hat)) / (min(t_end, t_hat + K) - t_hat)}
#'
#' \strong{7) AUC_deviation}
#'
#' Trapezoidal area under the curve of the milk deficit \code{baseline(t) - obs(t)}
#' across the whole episode. It summarizes how much milk was lost and for how long.
#'
#' Conceptually:
#'
#' \code{AUC_deviation = integral_{t_start..t_end} [baseline(t) - obs(t)] dt}
#'
#' In practice this is approximated via the trapezoidal rule on discrete DIMs.
#'
#' \strong{8) prod_decline_slope_amp}
#'
#' Product of the decline slope (anchor -> nadir) and the amplitude (anchor - nadir).
#' It combines speed and depth of the decline into a single indicator of how
#' "aggressive" the drop is.
#'
#' \code{prod_decline_slope_amp = DeclineSlope * A}
#'
#' \strong{9) prod_recovery_slope_TTB}
#'
#' Product of the recovery slope (nadir -> recovery) and time-to-baseline (TTB).
#' It combines how fast the animal recovers with how long recovery takes,
#' summarizing recovery efficiency.
#'
#' \code{prod_recovery_slope_TTB = RecoverySlope * TTB}
#'
#' @export
milkloss_detect <- function(
    data,
    id_col,
    dim_col,
    MY_col  = "MY_real",
    MY_pred,
    dim_start = 1L,
    dim_end   = 305L,
    rec_mode  = c("pctbase","band","resid"),
    drop_pct  = 0.10,
    min_len   = 1L,
    tol       = 0.05,
    stick     = 3L,
    rec       = 1.0
) {
  rec_mode <- match.arg(rec_mode)

  col_get <- function(d, nm) d[[nm]]
  slope_ols <- function(x, y) {
    x <- as.numeric(x); y <- as.numeric(y)
    if (length(x) < 2L) return(NA_real_)
    xv <- x - mean(x, na.rm = TRUE)
    denom <- sum(xv^2, na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) return(NA_real_)
    num <- sum(xv * (y - mean(y, na.rm = TRUE)), na.rm = TRUE)
    as.numeric(num / denom)
  }

  data <- data[order(col_get(data, id_col), col_get(data, dim_col)), , drop = FALSE]
  ids  <- unique(col_get(data, id_col))

  episodes_all <- list()
  aggs_all     <- list()

  for (cid in ids) {
    g <- data[col_get(data, id_col) == cid, , drop = FALSE]
    g <- g[stats::complete.cases(g[, c(dim_col, MY_col, MY_pred), drop = FALSE]), , drop = FALSE]
    g <- g[(col_get(g, dim_col) >= dim_start) & (col_get(g, dim_col) <= dim_end), , drop = FALSE]
    if (nrow(g) == 0L) next

    g <- g[order(col_get(g, dim_col)), , drop = FALSE]

    dims <- as.numeric(col_get(g, dim_col))
    obs  <- as.numeric(col_get(g, MY_col))
    base <- as.numeric(col_get(g, MY_pred))
    r    <- obs - base
    Tn   <- length(obs)

    ep_rows <- list()

    i <- 1L
    while (i <= Tn) {
      if (r[i] >= 0) { i <- i + 1L; next }

      p <- i - 1L
      while (p >= 1L && r[p] < 0) p <- p - 1L
      anchor_idx <- if (p >= 1L) p else NA_integer_
      anchor_ref <- if (!is.na(anchor_idx)) base[anchor_idx] else base[i]

      start_idx <- i
      neg_count <- 0L
      peak_drop_pct <- -Inf
      nadir_idx <- i
      nadir_r   <- r[i]

      recovered_row <- NA_integer_
      band_keep <- 0L; resid_keep <- 0L; pct_keep <- 0L

      j <- i
      while (j <= Tn && r[j] < 0) {
        neg_count <- neg_count + 1L

        if (r[j] < nadir_r) { nadir_r <- r[j]; nadir_idx <- j }

        drop_pct_j <- (anchor_ref - obs[j]) / max(anchor_ref, .Machine$double.eps)
        if (drop_pct_j > peak_drop_pct) peak_drop_pct <- drop_pct_j

        if (j > start_idx) {
          if (rec_mode == "band") {
            if (obs[j] >= base[j]*(1 - tol) && obs[j] <= base[j]*(1 + tol)) {
              band_keep <- band_keep + 1L
              if (band_keep >= max(1L, as.integer(stick))) {
                recovered_row <- j; break
              }
            } else band_keep <- 0L
          } else if (rec_mode == "resid") {
            thr <- rec * abs(nadir_r)
            if ((r[j] - nadir_r) >= thr) {
              resid_keep <- resid_keep + 1L
              if (resid_keep >= max(1L, as.integer(stick))) {
                recovered_row <- j; break
              }
            } else resid_keep <- 0L
          } else if (rec_mode == "pctbase") {
            if (obs[j] >= base[j] * rec) {
              pct_keep <- pct_keep + 1L
              if (pct_keep >= max(1L, as.integer(stick))) {
                recovered_row <- j; break
              }
            } else pct_keep <- 0L
          }
        }
        j <- j + 1L
      }

      if (is.na(recovered_row) && j <= Tn && r[j] >= 0) recovered_row <- j

      end_idx <- min(j, Tn)
      valid <- (neg_count >= as.integer(min_len)) && (peak_drop_pct >= drop_pct)
      if (!valid) {
        i <- if (!is.na(recovered_row)) recovered_row + 1L else max(end_idx + 1L, i + 1L)
        next
      }

      amp_anchor <- anchor_ref - obs[nadir_idx]

      aud <- 0
      for (t in start_idx:end_idx) {
        next_dim   <- if (t + 1L <= end_idx) dims[t + 1L] else dims[t] + 1
        delta_days <- max(1, as.integer(next_dim - dims[t]))
        aud <- aud + (base[t] - obs[t]) * delta_days
      }

      auc_dev <- 0
      if (start_idx == end_idx) {
        auc_dev <- base[start_idx] - obs[start_idx]
      } else {
        for (t in start_idx:(end_idx - 1L)) {
          x0 <- dims[t]
          x1 <- dims[t + 1L]
          delta_days <- max(1, as.integer(x1 - x0))
          y0 <- max(0, base[t] - obs[t])
          y1 <- max(0, base[t + 1L] - obs[t + 1L])
          auc_dev <- auc_dev + 0.5 * (y0 + y1) * delta_days
        }
      }

      days_to_nadir <- as.integer(dims[nadir_idx] - dims[start_idx])
      days_from_nadir_to_rec <- if (!is.na(recovered_row)) {
        as.integer(dims[recovered_row] - dims[nadir_idx])
      } else NA_integer_

      # TTB
      TTB_days <- days_from_nadir_to_rec

      # Half-life relative to anchor
      thalf_days <- NA_real_
      if (is.finite(amp_anchor) && amp_anchor > 1e-8) {
        half_level <- anchor_ref - amp_anchor / 2
        hit <- NA_integer_
        if (nadir_idx + 1L <= end_idx) {
          for (tt in (nadir_idx + 1L):end_idx) {
            if (obs[tt] >= half_level) { hit <- tt; break }
          }
        }
        if (!is.na(hit)) thalf_days <- as.numeric(dims[hit] - dims[nadir_idx])
      }

      decline_slope <- if (!is.na(anchor_idx) && dims[nadir_idx] != dims[anchor_idx]) {
        slope_ols(dims[anchor_idx:nadir_idx], obs[anchor_idx:nadir_idx])
      } else NA_real_

      recovery_slope <- if (!is.na(recovered_row) && dims[recovered_row] != dims[nadir_idx]) {
        slope_ols(dims[nadir_idx:recovered_row], obs[nadir_idx:recovered_row])
      } else NA_real_

      prod_decline_slope_amp <- decline_slope * amp_anchor
      prod_recovery_slope_TTB <- recovery_slope * TTB_days

      ep_rows[[length(ep_rows) + 1L]] <- data.frame(
        ID = cid,
        episode_index = length(ep_rows) + 1L,
        start_DIM = as.integer(dims[start_idx]),
        end_DIM   = as.integer(dims[end_idx]),
        nadir_DIM = as.integer(dims[nadir_idx]),
        recovery_DIM = if (!is.na(recovered_row)) as.integer(dims[recovered_row]) else NA_integer_,
        peak_drop_pct = as.numeric(peak_drop_pct),
        anchor_ref = as.numeric(anchor_ref),
        amplitude_anchor = as.numeric(amp_anchor),
        ML_per_event = as.numeric(aud),
        AUC_deviation = as.numeric(auc_dev),
        TTB_days = as.numeric(TTB_days),
        t_half_anchor_days = as.numeric(thalf_days),
        decline_slope = as.numeric(decline_slope),
        recovery_slope = as.numeric(recovery_slope),
        prod_decline_slope_amp = as.numeric(prod_decline_slope_amp),
        prod_recovery_slope_TTB = as.numeric(prod_recovery_slope_TTB),
        days_to_nadir = as.numeric(days_to_nadir),
        days_from_nadir_to_recovery = as.numeric(TTB_days),
        stringsAsFactors = FALSE
      )

      i <- if (!is.na(recovered_row)) recovered_row + 1L else max(end_idx + 1L, i + 1L)
    }

    episodes_data <- if (length(ep_rows)) {
      do.call(rbind, ep_rows)
    } else {
      data.frame(
        ID = character(0), episode_index = integer(0),
        start_DIM = integer(0), end_DIM = integer(0), nadir_DIM = integer(0),
        recovery_DIM = integer(0), peak_drop_pct = numeric(0),
        anchor_ref = numeric(0), amplitude_anchor = numeric(0),
        ML_per_event = numeric(0), AUC_deviation = numeric(0),
        TTB_days = numeric(0),
        t_half_anchor_days = numeric(0),
        decline_slope = numeric(0),
        recovery_slope = numeric(0),
        prod_decline_slope_amp = numeric(0),
        prod_recovery_slope_TTB = numeric(0),
        days_to_nadir = numeric(0), days_from_nadir_to_recovery = numeric(0),
        stringsAsFactors = FALSE
      )
    }

    episodes_all[[length(episodes_all) + 1L]] <- episodes_data

    aggs_all[[length(aggs_all) + 1L]] <- data.frame(
      ID = cid,
      n_episodes = nrow(episodes_data),
      total_ML   = if (nrow(episodes_data)) sum(episodes_data$ML_per_event,   na.rm = TRUE) else 0,
      total_auc_dev = if (nrow(episodes_data)) sum(episodes_data$AUC_deviation, na.rm = TRUE) else 0,
      total_days = if (nrow(episodes_data)) sum(episodes_data$end_DIM - episodes_data$start_DIM, na.rm = TRUE) else 0,
      mean_amp   = if (nrow(episodes_data)) mean(episodes_data$amplitude_anchor, na.rm = TRUE) else NA_real_,
      mean_thalf = if (nrow(episodes_data)) mean(episodes_data$t_half_anchor_days, na.rm = TRUE) else NA_real_,
      median_ttb = if (nrow(episodes_data)) stats::median(episodes_data$TTB_days, na.rm = TRUE) else NA_real_,
      mean_prod_decline_slope_amp =
        if (nrow(episodes_data)) mean(episodes_data$prod_decline_slope_amp, na.rm = TRUE) else NA_real_,
      mean_prod_recovery_slope_TTB =
        if (nrow(episodes_data)) mean(episodes_data$prod_recovery_slope_TTB, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  episodes_out <- if (length(episodes_all)) {
    do.call(rbind, episodes_all)
  } else {
    data.frame(
      ID = character(0), episode_index = integer(0),
      start_DIM = integer(0), end_DIM = integer(0), nadir_DIM = integer(0),
      recovery_DIM = integer(0), peak_drop_pct = numeric(0),
      anchor_ref = numeric(0), amplitude_anchor = numeric(0),
      ML_per_event = numeric(0), AUC_deviation = numeric(0),
      TTB_days = numeric(0),
      t_half_anchor_days = numeric(0),
      decline_slope = numeric(0),
      recovery_slope = numeric(0),
      prod_decline_slope_amp = numeric(0),
      prod_recovery_slope_TTB = numeric(0),
      days_to_nadir = numeric(0), days_from_nadir_to_recovery = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  aggregates_out <- if (length(aggs_all)) {
    do.call(rbind, aggs_all)
  } else {
    data.frame(
      ID = character(0), n_episodes = integer(0),
      total_ML = numeric(0), total_auc_dev = numeric(0),
      total_days = numeric(0),
      mean_amp = numeric(0), mean_thalf = numeric(0),
      median_ttb = numeric(0),
      mean_prod_decline_slope_amp = numeric(0),
      mean_prod_recovery_slope_TTB = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  list(episodes = episodes_out, aggregates = aggregates_out)
}
