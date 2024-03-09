#' Helper for plotting
#'
#' @param Gibbs Gibbs output object
#' @param burnin `integer` giving number of burn in periods
#' @param path character giving the path to the where plots are stored
#' @param nameMat names of matrices
#' @param nameReg names of regressors
#' @param njointfac number of joint factors
#' @param type e.g. "allidio"
#' @param predictionCI logical; if `TRUE`, the prediction intervalls are added
#' @param onlyY some setting
#' @param hierachPrior if there are hierarchical prior
#'
#' @return side effect function returning the diagnostic plots
#' @export
save_Gibbs_plots <- function(Gibbs, burnin, path, nameMat, nameReg, njointfac,
                             type, predictionCI = TRUE, onlyY = FALSE,
                             hierachPrior, true_values = NULL) {
  if (is.null(Gibbs$errorMsg)) {
    yObs <- Gibbs$initials$yObs
    npara <- Gibbs$initials$npara
    nreg <- Gibbs$initials$nreg
    TT <- dim(yObs)[2]
    N <- dim(yObs)[1] / npara
    if (missing(nameMat)) {
      nameMat <- NULL
      countries <- 1:N
    } else {
      countries <- unique(nameMat[, 1])
    }
    if (missing(nameReg)) {
      nameReg <- NULL
    }
    sampleA <- Gibbs$initials$sampleA
    countryA <- Gibbs$initials$countryA
    sampleV <- Gibbs$initials$sampleV
    diagA <- Gibbs$initials$diagA
    plotlist <- plot_Gibbs_seq(Gibbs, burnin, nameMat, nameReg, njointfac, type, predictionCI = predictionCI, hierachPrior = hierachPrior)
    B0 <- round(Gibbs$initials$B0[1, 1, 1], 2)
    D0 <- Gibbs$initials$D0[1, 1, 1]
    Omega0 <- Gibbs$initials$Omega0[1, 1]
    if (nreg != 0) {
      dimOmega0 <- dim(Gibbs$initials$Omega0)[1]
      OmegaD0 <- Gibbs$initials$Omega0[dimOmega0, dimOmega0]
    } else {
      OmegaD0 <- NULL
    }
    incObsNew <- Gibbs$initials$incObsNew
    p_joint <- njointfac
    itermax <- Gibbs$initials$itermax
    itermax <- itermax - burnin

    path <- paste0(path, "/IO", incObsNew, paste0(nameReg, collapse = "_"))
    dir.create(path)

    covScale <- Gibbs$initials$covScale

    B_plot <- plotlist$B_plot
    y_plot <- plotlist$y_plot
    f_plot <- plotlist$f_plot
    D_plot <- plotlist$D_plot
    f_Gibbs_seq_plot <- plotlist$f_Gibbs_seq_plot
    Bf_Gibbs_seq_plot <- plotlist$Bf_Gibbs_seq_plot
    Mu_plot <- plotlist$Mu_plot
    Omega_plot <- plotlist$Omega_plot


    if (hierachPrior) {
      mub <- Gibbs$initials$mu_b0[1]
      Sigmab <- Gibbs$initials$Sigma_b0[1, 1]
      alphab <- Gibbs$initials$alpha_b0[1]
      betab <- Gibbs$initials$beta_b0[1]
      path_file_hierach <- paste0("_alphb", alphab, "_betb", betab, "_Sigb", Sigmab, "_mub", mub)
    } else {
      path_file_hierach <- NULL
    }


    if (sampleV) {
      V_plot <- plotlist$V_plot

      alpha0 <- Gibbs$initials$alpha0
      beta0 <- Gibbs$initials$beta0
      path_adj <- paste0(path, "/pj", njointfac, "_B", B0, "_Om", Omega0, "_D", D0, "_OmD", OmegaD0, "_alpha", alpha0, "_beta", beta0, "_IO", incObsNew, path_file_hierach)
      dir.create(path_adj)
      file_name <- paste0("_pj", njointfac, "_B", B0, "_Om", Omega0, "_D", D0, "_OmD", OmegaD0, "_IO", incObsNew, path_file_hierach)
      plot_top <- paste(
        itermax, "draws", "(burnin =", burnin, "),", "T =", TT, ", N =", N, ", common factors =", p_joint, ", B0 =", B0, ", Omega0 =", Omega0, ", D0 =", D0, ", OmegaD0 =", OmegaD0,
        ", alpha0 =", alpha0, ", beta0 =", beta0
      )
      if (!onlyY) {
        save_zoom(plot_list = V_plot[1:floor(length(V_plot) / 2)], as.table = F, ncol = N * npara / 6, path = path_adj, file = paste0("V_part1", file_name), top = paste("Variances (diagonal)", plot_top))
        save_zoom(plot_list = V_plot[floor(length(V_plot) / 2 + 1):length(V_plot)], as.table = F, ncol = N * npara / 6, path = path_adj, file = paste0("V_part2", file_name), top = paste("Variances (diagonal)", plot_top))
      }
    } else if (sampleA) {
      startA <- Gibbs$initials$A

      startA_lt <- c(startA[lower.tri(startA, diag = T)])
      A_plot <- plotlist$A_plot
      A_name <- startA[1, 1]
      if (diagA) {
        shape0 <- Gibbs$initials$shape0
        rate0 <- Gibbs$initials$rate0
        path_adj <- paste0(path, "/B", B0, "_Om", Omega0, "_D", D0, "_OmD", OmegaD0, "_A", A_name, "_shape", shape0, "_rate", rate0, "_IO", incObsNew, path_file_hierach)
        file_name <- paste0("_B", B0, "_Om", Omega0, "_D", D0, "_OmD", OmegaD0, "_A", A_name, "_shape", shape0, "_rate", rate0, "_IO", incObsNew, path_file_hierach)
        plot_top <- paste(itermax, "draws", "(burnin =", burnin, "),", "T =", TT, ", N =", N, ", common factors =", p_joint, ", B0 =", B0, ", Omega0 =", Omega0, ", D0 =", D0, ", OmegaD0 =", OmegaD0, ", Shape =", shape0, ", Rate =", rate0)
      } else {
        Psi0 <- Gibbs$initials$Psi0[1, 1]
        nu0 <- Gibbs$initials$nu0
        path_adj <- paste0(path, "/pj", njointfac, "_cs", covScale, "_B", B0, "_Om", Omega0, "_D", D0, "_OmD", OmegaD0, "_A", A_name, "_Psi", Psi0, "_nu", nu0, "_IO", incObsNew, path_file_hierach)
        file_name <- paste0("_pj", njointfac, "_cs", covScale, "_B", B0, "_Om", Omega0, "_D", D0, "_OmD", OmegaD0, "_A", A_name, "_Psi", Psi0, "_nu", nu0, "_IO", incObsNew, path_file_hierach)
        plot_top <- paste(itermax, "draws", "(burnin =", burnin, "),", "T =", TT, ", N =", N, ", common factors =", p_joint, ", B0 =", B0, ", Omega0 =", Omega0, ", D0 =", D0, ", OmegaD0 =", OmegaD0, ", Psi0 =", Psi0, ", nu0 =", nu0)
      }
      dir.create(path_adj)

      if (countryA) {
        Afile_name <- lapply(countries, \(x) paste0("_", x, "_pj", njointfac, "_cs", covScale, "_B", B0, "_Om", Omega0, "_D", D0, "_OmD", OmegaD0, "_A", A_name, "_Psi", Psi0, "_nu", nu0, "_IO", incObsNew))
        Aplot_top <- lapply(countries, \(x) paste(x, ",", itermax, "draws", "(burnin =", burnin, "),", "T =", TT, ", N =", N, ", common factors =", p_joint, ", B0 =", B0, ", Omega0 =", Omega0, ", D0 =", D0, ", OmegaD0 =", OmegaD0, ", Psi0 =", Psi0, ", nu0 =", nu0))
        if (!onlyY) {
          for (i in 1:N) {
            plot_num <- npara * (npara + 1) / 2
            save_zoom(plot_list = A_plot[((i - 1) * plot_num + 1):(i * plot_num)], as.table = F, ncol = 2, path = path_adj, file = paste0("A_lowerdiag", Afile_name[[i]][1]), top = paste("Adjustment matrix (lower diagonal)", Aplot_top[[i]][1]))
          }
        }
      } else {
        if (!onlyY) {
          save_zoom(plot_list = A_plot, as.table = F, ncol = 2, path = path_adj, file = paste0("A_lowerdiag", file_name), top = paste("Adjustment matrix (lower diagonal)", plot_top))
        }
      }
    } else {
      path_adj <- paste0(path, "/pj", njointfac, "_cs", covScale, "_B", B0, "_Om", Omega0, "_D", D0, "_OmD", OmegaD0, "_IO", incObsNew, path_file_hierach)
      dir.create(path_adj)
      file_name <- paste0("_pj", njointfac, "_cs", covScale, "_B", B0, "_Om", Omega0, "_D", D0, "_OmD", OmegaD0, "_IO", incObsNew, path_file_hierach)
      plot_top <- paste(itermax, "draws", "(burnin =", burnin, "),", "T =", TT, ", N =", N, ", common factors =", p_joint, ", cs =", covScale, ", B0 =", B0, ", Omega0 =", Omega0, ", D0 =", D0, ", OmegaD0 =", OmegaD0)
    }




    if (!onlyY) {
      if (njointfac != 0) {
        nfiles_joint <- floor(2 * N * npara * njointfac / 30)
        for (i in 1:nfiles_joint) {
          save_zoom(B_plot[(floor((2 * N * npara * njointfac) / nfiles_joint) * (i - 1) + 1):(floor((2 * N * npara * njointfac) / nfiles_joint) * i)],
            ncol = 5, as.table = F, path = path_adj,
            file = paste0("Loadings_common_part", i, file_name), top = paste("Loadings (common),", plot_top)
          )
        }

        if (nfiles_joint != (2 * N * npara * njointfac / 30)) {
          save_zoom(B_plot[(floor((2 * N * npara * njointfac) / nfiles_joint) * nfiles_joint + 1):(2 * N * npara * njointfac)],
            ncol = 5, as.table = F, path = path_adj,
            file = paste0("Loadings_common_part", nfiles_joint + 1, file_name), top = paste("Loadings (common),", plot_top)
          )
        }
      } else {
        nfiles_joint <- 0
      }

      nfiles <- floor((length(B_plot) - 2 * N * npara * njointfac) / 30)

      for (i in 1:nfiles) {
        save_zoom(B_plot[(floor(length(B_plot) / (nfiles + nfiles_joint)) * (i - 1) + 1 + 2 * N * npara * njointfac):(floor(length(B_plot) / (nfiles + nfiles_joint)) * i + 2 * N * npara * njointfac)],
          ncol = 5, as.table = F, path = path_adj,
          file = paste0("Loadings_idio_part", i, file_name), top = paste("Loadings (idiosyncratic),", plot_top)
        )
      }

      if (nfiles != (length(B_plot) - 2 * N * npara * njointfac) / 30) {
        save_zoom(B_plot[(floor(length(B_plot) / (nfiles + nfiles_joint)) * nfiles + 1 + 2 * N * npara * njointfac):length(B_plot)],
          ncol = 5, as.table = F, path = path_adj,
          file = paste0("Loadings_idio_part", nfiles + 1, file_name), top = paste("Loadings (idiosyncratic),", plot_top)
        )
      }

      if (nreg != 0) {
        nfiles <- floor(length(D_plot) / 30)

        for (i in 1:nfiles) {
          save_zoom(D_plot[(floor(length(D_plot) / nfiles) * (i - 1) + 1):(floor(length(D_plot) / nfiles) * i)],
            ncol = 5, as.table = F, path = path_adj,
            file = paste0("Macros_part", i, file_name), top = paste("Macros (partial effects),", plot_top)
          )
        }

        if (nfiles != length(D_plot) / 30) {
          save_zoom(D_plot[(floor(length(D_plot) / nfiles) * nfiles + 1):length(D_plot)],
            ncol = 5, as.table = F, path = path_adj,
            file = paste0("Macros_part", nfiles + 1, file_name), top = paste("Macros (partial effects),", plot_top)
          )
        }
      }


      save_zoom(f_plot, ncol = N * npara / 5, path = path_adj, file = paste0("f", file_name), top = paste("Factors (mean of Gibbs draws),", plot_top))
      save_zoom(f_Gibbs_seq_plot, ncol = 4, path = path_adj, file = paste0("f_Gibbs_seq", file_name), top = paste("Factors Gibbs seq.,", plot_top))
      save_zoom(Bf_Gibbs_seq_plot, ncol = 4, path = path_adj, file = paste0("Bf_Gibbs_seq", file_name), top = paste("Bf Gibbs seq.,", plot_top))

      if (hierachPrior) {
        save_zoom(Mu_plot, ncol = 4, as.table = F, path = path_adj, file = paste0("Mu", file_name), top = paste("Prior Mean Gibbs seq.", plot_top))
        save_zoom(Omega_plot, ncol = 4, as.table = F, path = path_adj, file = paste0("Omega", file_name), top = paste("Prior Variance Gibbs seq.", plot_top))
      }
    }
    save_zoom(y_plot, ncol = N * npara / 5, path = path_adj, file = paste0("y", "_CI", substr(as.character(predictionCI), 1, 1), file_name), top = paste("Data and Prediction,", plot_top))
  }
}
plot_Gibbs_seq <- function(Gibbs,
                           burnin,
                           nameMat,
                           nameReg,
                           njointfac,
                           type,
                           predictionCI,
                           onlyY = FALSE,
                           hierachPrior) {
  npara <- Gibbs$initials$npara
  nreg <- Gibbs$initials$nreg
  sampleA <- Gibbs$initials$sampleA
  sampleV <- Gibbs$initials$sampleV
  countryA <- Gibbs$initials$countryA

  yObs <- Gibbs$initials$yObs
  TT <- dim(yObs)[2]
  N <- dim(yObs)[1] / npara
  p_joint <- njointfac
  p <- dim(Gibbs$B)[2]
  itermax <- Gibbs$initials$itermax

  itermax <- itermax - burnin
  selecetedObs <- (burnin + 1):(itermax + burnin)

  if (nreg != 0) {
    D_Gibbs_mean <- apply(Gibbs$D[, , selecetedObs], c(1, 2), mean)
    wReg <- Gibbs$initials$wReg
    regFit <- t(D_Gibbs_mean %*% wReg)
    yfit_mat <- sapply(1:dim(Gibbs$B)[3], \(x) Gibbs$B[, , x] %*% Gibbs$f[, , x] + Gibbs$D[, , x] %*% wReg)
  } else {
    regFit <- matrix(0, ncol = N * npara, nrow = TT)
    yfit_mat <- sapply(1:dim(Gibbs$B)[3], \(x) Gibbs$B[, , x] %*% Gibbs$f[, , x])
  }

  yfit_array <- array(yfit_mat, dim = c(N * npara, TT, itermax + burnin))

  f_Gibbs_mean <- apply(Gibbs$f[, , selecetedObs], c(1, 2), mean)
  B_Gibbs_mean <- apply(Gibbs$B[, , selecetedObs], c(1, 2), mean)
  yfit_mean <- apply(yfit_array[, , selecetedObs], c(1, 2), mean)
  if (predictionCI) {
    yfit_uci <- apply(yfit_array[, , selecetedObs], c(1, 2), stats::quantile, probs = 0.975)
    yfit_lci <- apply(yfit_array[, , selecetedObs], c(1, 2), stats::quantile, probs = 0.025)
  } else {
    yfit_uci <- NULL
    yfit_lci <- NULL
  }

  y_plot <- vector("list", npara * N)
  for (i in 1:(npara * N)) {
    y_Gibbs <- dplyr::tibble(Time = 1:TT, y_data = yObs[i, ], y_Gibbs = yfit_mean[i, ], y_Gibbs_uci = yfit_uci[i, ], y_Gibbs_lci = yfit_lci[i, ], y_Gibbs_Bf = t(f_Gibbs_mean) %*% t(B_Gibbs_mean)[, i] + regFit[, i])
    # y_Gibbs <- dplyr::tibble(Time = 1:TT, y_data = yObs[i,], y_Gibbs =  t(f_Gibbs_mean)%*%t(B_Gibbs_mean)[,i] )
    y_plot[[i]] <- local({
      i <- i
      if (is.null(nameMat)) {
        ylab_text <- paste("y", "Obs", i)
      } else {
        ylab_text <- paste(nameMat[i, 1], nameMat[i, 2])
      }

      if (i == npara * N) {
        legend_ <- "bottom"
      } else {
        legend_ <- "none"
      }
      ggplot2::ggplot(data = y_Gibbs, mapping = ggplot2::aes(x = Time, y = y_data, color = "y", linetype = "line")) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::geom_line(ggplot2::aes(x = Time, y = y_Gibbs, color = "y fit", linetype = "line")) +
        ggplot2::geom_line(ggplot2::aes(x = Time, y = y_Gibbs_Bf, color = "y Bayes Pred.", linetype = "line")) +
        {
          if (predictionCI) ggplot2::geom_line(ggplot2::aes(x = Time, y = y_Gibbs_uci, color = "y fit", linetype = "cf_band"))
        } +
        {
          if (predictionCI) ggplot2::geom_line(ggplot2::aes(x = Time, y = y_Gibbs_lci, color = "y fit", linetype = "cf_band"))
        } +
        {
          if (predictionCI) ggplot2::scale_linetype_manual(values = c("dashed", "solid"))
        } +
        {
          if (!predictionCI) ggplot2::scale_linetype_manual(values = c("solid"))
        } +
        ggplot2::guides(linetype = "none") +
        ggplot2::ylab(ylab_text) +
        ggplot2::labs(color = "", title = "") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = legend_, axis.text.x = ggplot2::element_text(size = 8)) +
        ggplot2::scale_x_continuous(breaks = 1:TT)
    })
  }

  if (!onlyY) {
    Bf_Gibbs_seq_i <- 1:4
    Bf_Gibbs_seq_t <- c(1, seq(floor(TT / 3), floor(TT / 3) * 3, floor(TT / 3)))
    Bf_Gibbs_seq_plot <- vector("list", length(Bf_Gibbs_seq_i) * length(Bf_Gibbs_seq_t))
    j <- 1
    for (s in Bf_Gibbs_seq_i) {
      for (t in Bf_Gibbs_seq_t) {
        Bf_Gibbs <- dplyr::tibble(Gibbs_draws = 1:length(selecetedObs), Gibbs_Bf = yfit_array[s, t, selecetedObs])

        Bf_Gibbs_seq_plot[[j]] <- local({
          j <- j
          s <- s
          t <- t
          # if(p_joint == 0){s2 <- s2 + 1}
          # if(s %in% 1:p_joint & p_joint != 0){
          #   title_y <- paste0("","p = J",s,", T =",t)
          # }else{
          #   title_y <- paste0("","p = I",s2-1,", T =",t)
          # }
          if (is.null(nameMat)) {
            ylab_text <- paste0("i = ", s, ", t =", t)
          } else {
            ylab_text <- paste(nameMat[s, 1], nameMat[s, 2], "t =", t)
          }
          ggplot2::ggplot(data = Bf_Gibbs, ggplot2::aes(x = Gibbs_draws, y = Gibbs_Bf, color = "Gibbs draws")) +
            ggplot2::geom_line() +
            ggplot2::scale_color_manual(values = c("darkgrey")) +
            ggplot2::ylab(ylab_text) +
            ggplot2::labs(color = "", title = paste("Mean(Gibbs) =", round(mean(Bf_Gibbs$Gibbs_Bf), 2), ", acf1 =", round(acf(Bf_Gibbs$Gibbs_Bf, plot = F)[[1]][2], 2))) +
            ggplot2::theme_bw()
        })
        j <- j + 1
      }
    }
    # f_Gibbs_mean <- apply(Gibbs$f[,,selecetedObs],c(1,2),mean)
    f_plot <- vector("list", p)
    j <- 1
    for (i in c(2:p, 1)) {
      f_Gibbs <- dplyr::tibble(Time = 1:TT, Gibbs_mean_f = t(f_Gibbs_mean)[, i])
      # Gibbs_mean_f = Gibbs$f[i,])

      f_plot[[j]] <- local({
        j <- j
        i <- i
        if (p_joint != 0) {
          if (i == 1) {
            title_ <- "Common"
          } else {
            if (is.null(nameMat)) {
              title_ <- paste("Idio", i - 1)
            } else {
              title_ <- paste("Idio", nameMat[j, 1], nameMat[j, 2])
            }
          }
        } else {
          title_ <- paste(i)
        }
        if (i == p) {
          legend_ <- "bottom"
        } else {
          legend_ <- "none"
        }
        ggplot2::ggplot(data = f_Gibbs, ggplot2::aes(x = Time, y = Gibbs_mean_f, color = paste0("Mean of Gibbs draws (", paste(itermax), ")"))) +
          ggplot2::geom_line() +
          ggplot2::ylab(paste(title_)) +
          ggplot2::labs(color = "", title = "") +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position = "none")
        # scale_x_continuous(breaks = 1:TT)+
        # theme(legend.position = "none", axis.text.x = element_text(size = 4.5))
      })
      j <- j + 1
    }


    if (type == "allidio") {
      selectR_complete_vec <- c(rep(1, N * npara * p_joint + 1), rep(c(rep(0, N * npara), 1), N * npara - 1))
      selectR_complete <- diag(selectR_complete_vec)[-which(selectR_complete_vec == 0), ]
      # selectC <- rep(0, npara + npara * p_joint)
    } else if (type == "countryidio") {
      selectR_complete_vec <- c(rep(1, N * npara + npara), rep(c(rep(0, N * npara), rep(1, npara)), N - 1))
      selectR_complete <- diag(selectR_complete_vec)[-which(selectR_complete_vec == 0), ]
    } else if (type == "countryidio_nomu") {
      selectR_complete_vec <- c(rep(1, N * npara + npara - 1), rep(c(rep(0, N * npara + 1), rep(1, npara - 1)), N - 1), 0)
      selectR_complete <- diag(selectR_complete_vec)[-which(selectR_complete_vec == 0), ]
    }

    Gibbs_B_nonzero <- apply(Gibbs$B[, , selecetedObs], 3, function(x) {
      selectR_complete %*% c(x)
    })
    GibbsB <- Gibbs_B_nonzero
    GibbsB_mean <- apply(GibbsB, 1, mean)
    GibbsB_uci <- apply(GibbsB, 1, stats::quantile, probs = 0.975)
    GibbsB_lci <- apply(GibbsB, 1, stats::quantile, probs = 0.025)


    # B_plot <- vector("list", N * (npara + npara * p_joint))
    B_plot <- vector("list", 2 * sum(selectR_complete_vec))
    nameMat_ext <- rbind(nameMat, nameMat)
    # for(i in 1:(N * (npara + npara * p_joint))){
    j <- 1
    for (i in 1:sum(selectR_complete_vec)) {
      B_Gibbs <- dplyr::tibble(Gibbs_draws = 1:length(GibbsB[1, ]), Gibbs_c = GibbsB[i, ])
      B_plot[[j]] <-
        local({
          i <- i
          if (is.null(nameMat)) {
            ylab_text <- paste("Load", i)
          } else {
            ylab_text <- paste(nameMat_ext[i, 1], nameMat_ext[i, 2])
          }
          BayEstB <- GibbsB_mean[i]
          uci <- round(GibbsB_uci[i], 2)
          lci <- round(GibbsB_lci[i], 2)
          ggplot2::ggplot(data = B_Gibbs, mapping = ggplot2::aes(x = Gibbs_draws, y = Gibbs_c, color = "Gibbs draws", linetype = "line")) +
            ggplot2::geom_line() +
            ggplot2::geom_hline(ggplot2::aes(yintercept = BayEstB, color = "Bay. Est.", linetype = "line")) +
            ggplot2::geom_hline(ggplot2::aes(yintercept = lci, color = "Bay. Est.", linetype = "cf_band")) +
            ggplot2::geom_hline(ggplot2::aes(yintercept = uci, color = "Bay. Est.", linetype = "cf_band")) +
            ggplot2::scale_linetype_manual(values = c("dashed", "solid")) +
            # ggplot2::labs(color="", title = paste("Est =",round(BayEstB,3),", acf1 =", round(acf(B_Gibbs$Gibbs_c,plot=F)[[1]][2],3) ))+
            ggplot2::labs(color = "", title = paste0("Est =", round(BayEstB, 3), ", 95%-CI = ", "[", lci, ", ", uci, "]")) +
            ggplot2::ylab(ylab_text) +
            ggplot2::xlab("Gibbs draws") +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position = "none", text = ggplot2::element_text(size = 12))
        })

      j <- j + 1
      B_plot[[j]] <- local({
        i <- i
        if (is.null(nameMat)) {
          ylab_text <- paste("ACF Load", i)
        } else {
          ylab_text <- paste("ACF", nameMat_ext[i, 1], nameMat_ext[i, 2])
        }
        bacf <- acf(B_Gibbs$Gibbs_c, plot = FALSE)
        bacfdf <- with(bacf, data.frame(lag, acf))
        ggplot2::ggplot(data = bacfdf, mapping = ggplot2::aes(x = lag, y = acf)) +
          ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
          ggplot2::geom_segment(mapping = ggplot2::aes(xend = lag, yend = 0)) +
          ggplot2::ylab(ylab_text)
      })
      j <- j + 1
    }


    if (nreg != 0) {
      Dones_vec <- c(makeDstart(npara = npara, N = N, nreg = nreg, D_par = 1)$Dstack)
      selectR_D_complete <- diag(Dones_vec)[-which(Dones_vec == 0), ]
      Gibbs_D_nonzero <- apply(Gibbs$D[, , selecetedObs], 3, function(x) {
        selectR_D_complete %*% c(x)
      })
      GibbsD <- Gibbs_D_nonzero
      GibbsD_mean <- apply(GibbsD, 1, mean)
      GibbsD_uci <- apply(GibbsD, 1, stats::quantile, probs = 0.975)
      GibbsD_lci <- apply(GibbsD, 1, stats::quantile, probs = 0.025)

      D_plot <- vector("list", 2 * sum(selectR_D_complete))
      if (is.null(nameReg)) {
        nameReg <- as.character(1:nreg)
      }
      if (is.null(nameMat)) {
        nameMat <- as.character(1:(N * npara))
        nameMatReg <- matrix(rep(nameMat, each = nreg), ncol = 1)
      } else {
        nameMatReg <- matrix(rep(nameMat, each = nreg), ncol = 2)
      }
      nameMatReg_ext <- cbind(nameMatReg, nameReg)

      # for(i in 1:(N * (npara + npara * p_joint))){
      j <- 1
      for (i in 1:sum(selectR_D_complete)) {
        D_Gibbs <- dplyr::tibble(Gibbs_draws = 1:length(GibbsD[1, ]), Gibbs_c = GibbsD[i, ])
        D_plot[[j]] <-
          local({
            i <- i
            if (is.null(nameMat)) {
              ylab_text <- paste0(nameMatReg_ext[i, 1], nameMatReg_ext[i, 2])
            } else {
              ylab_text <- paste(nameMatReg_ext[i, 1], nameMatReg_ext[i, 2], nameMatReg_ext[i, 3])
            }
            BayEstD <- GibbsD_mean[i]
            uci <- round(GibbsD_uci[i], 2)
            lci <- round(GibbsD_lci[i], 2)
            ggplot2::ggplot(data = D_Gibbs, mapping = ggplot2::aes(x = Gibbs_draws, y = Gibbs_c, color = "Gibbs draws", linetype = "line")) +
              ggplot2::geom_line() +
              ggplot2::geom_hline(ggplot2::aes(yintercept = BayEstD, color = "Bay. Est.", linetype = "line")) +
              ggplot2::geom_hline(ggplot2::aes(yintercept = lci, color = "Bay. Est.", linetype = "cf_band")) +
              ggplot2::geom_hline(ggplot2::aes(yintercept = uci, color = "Bay. Est.", linetype = "cf_band")) +
              # ggplot2::labs(color="", title = paste("Est =",round(BayEstD, 3),", acf1 =", round(acf(D_Gibbs$Gibbs_c, plot=F)[[1]][2],3) ))+
              ggplot2::labs(color = "", title = paste0("Est =", round(BayEstD, 3), ", 95%-CI = ", "[", lci, ", ", uci, "]")) +
              ggplot2::scale_linetype_manual(values = c("dashed", "solid")) +
              ggplot2::ylab(ylab_text) +
              ggplot2::xlab("Gibbs draws") +
              ggplot2::theme_bw() +
              ggplot2::theme(legend.position = "none", text = ggplot2::element_text(size = 12))
          })


        j <- j + 1
        D_plot[[j]] <- local({
          i <- i
          if (is.null(nameMat)) {
            ylab_text <- paste("ACF Load", paste0(nameMatReg_ext[i, 1], nameMatReg_ext[i, 2]))
          } else {
            ylab_text <- paste("ACF", nameMatReg_ext[i, 1], nameMatReg_ext[i, 2], nameMatReg_ext[i, 3])
          }
          dacf <- acf(D_Gibbs$Gibbs_c, plot = FALSE)
          dacfdf <- with(dacf, data.frame(lag, acf))
          ggplot2::ggplot(data = dacfdf, mapping = ggplot2::aes(x = lag, y = acf)) +
            ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
            ggplot2::geom_segment(mapping = ggplot2::aes(xend = lag, yend = 0)) +
            ggplot2::ylab(ylab_text)
        })
        j <- j + 1
      }
    } else {
      D_plot <- NULL
    }


    if (hierachPrior) {
      GibbsMu <- Gibbs$BD0STORE
      GibbsOmega <- Gibbs$Omega0STORE

      n_plot <- nrow(GibbsMu)
      Mu_plot <- vector("list", n_plot)
      # for(i in 1:(N * (npara + npara * p_joint))){
      j <- 1
      for (i in 1:n_plot) {
        Mu_Gibbs <- dplyr::tibble(Gibbs_draws = 1:length(GibbsMu[1, ]), Gibbs_c = GibbsMu[i, ])
        Mu_plot[[j]] <-
          local({
            i <- i
            BayEstMu <- mean(Mu_Gibbs$Gibbs_c)
            ggplot2::ggplot(data = Mu_Gibbs, mapping = ggplot2::aes(x = Gibbs_draws, y = Gibbs_c, color = "Gibbs draws", linetype = "line")) +
              ggplot2::geom_line() +
              ggplot2::geom_hline(ggplot2::aes(yintercept = BayEstMu, color = "Bay. Est.", linetype = "line")) +
              # ggplot2::labs(color="", title = paste("Est =",round(BayEstB,3),", acf1 =", round(acf(B_Gibbs$Gibbs_c,plot=F)[[1]][2],3) ))+
              ggplot2::labs(color = "", title = paste0("Est =", round(BayEstMu, 3))) +
              ggplot2::ylab(paste(i)) +
              ggplot2::xlab("Gibbs draws") +
              ggplot2::theme_bw() +
              ggplot2::theme(legend.position = "none", text = ggplot2::element_text(size = 12))
          })

        j <- j + 1
        Mu_plot[[j]] <- local({
          i <- i
          if (is.null(nameMat)) {
            ylab_text <- paste("ACF", i)
          } else {
            ylab_text <- paste("ACF", i)
          }
          bacf <- acf(Mu_Gibbs$Gibbs_c, plot = FALSE)
          bacfdf <- with(bacf, data.frame(lag, acf))
          ggplot2::ggplot(data = bacfdf, mapping = ggplot2::aes(x = lag, y = acf)) +
            ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
            ggplot2::geom_segment(mapping = ggplot2::aes(xend = lag, yend = 0)) +
            ggplot2::ylab(ylab_text)
        })
        j <- j + 1
      }

      Omega_plot <- vector("list", n_plot)
      # for(i in 1:(N * (npara + npara * p_joint))){
      j <- 1
      for (i in 1:n_plot) {
        Omega_Gibbs <- dplyr::tibble(Gibbs_draws = 1:length(GibbsOmega[1, ]), Gibbs_c = GibbsOmega[i, ])
        Omega_plot[[j]] <-
          local({
            i <- i
            BayEstOmega <- mean(Omega_Gibbs$Gibbs_c)
            ggplot2::ggplot(data = Omega_Gibbs, mapping = ggplot2::aes(x = Gibbs_draws, y = Gibbs_c, color = "Gibbs draws", linetype = "line")) +
              ggplot2::geom_line() +
              ggplot2::geom_hline(ggplot2::aes(yintercept = BayEstOmega, color = "Bay. Est.", linetype = "line")) +
              # ggplot2::labs(color="", title = paste("Est =",round(BayEstB,3),", acf1 =", round(acf(B_Gibbs$Gibbs_c,plot=F)[[1]][2],3) ))+
              ggplot2::labs(color = "", title = paste0("Est =", round(BayEstOmega, 3))) +
              ggplot2::ylab(paste(i)) +
              ggplot2::xlab("Gibbs draws") +
              ggplot2::theme_bw() +
              ggplot2::theme(legend.position = "none", text = ggplot2::element_text(size = 12))
          })

        j <- j + 1
        Omega_plot[[j]] <- local({
          i <- i
          if (is.null(nameMat)) {
            ylab_text <- paste("ACF", i)
          } else {
            ylab_text <- paste("ACF", i)
          }
          bacf <- acf(Omega_Gibbs$Gibbs_c, plot = FALSE)
          bacfdf <- with(bacf, data.frame(lag, acf))
          ggplot2::ggplot(data = bacfdf, mapping = ggplot2::aes(x = lag, y = acf)) +
            ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
            ggplot2::geom_segment(mapping = ggplot2::aes(xend = lag, yend = 0)) +
            ggplot2::ylab(ylab_text)
        })
        j <- j + 1
      }
    } else {
      Mu_plot <- NULL
      Omega_plot <- NULL
    }

    f_Gibbs_seq_p <- 1:4
    f_Gibbs_seq_t <- c(1, seq(floor(TT / 3), floor(TT / 3) * 3, floor(TT / 3)))
    f_Gibbs_seq_plot <- vector("list", length(f_Gibbs_seq_p) * length(f_Gibbs_seq_t))
    j <- 1
    for (s in f_Gibbs_seq_p) {
      for (t in f_Gibbs_seq_t) {
        f_Gibbs <- dplyr::tibble(Gibbs_draws = 1:length(selecetedObs), Gibbs_f = Gibbs$f[s, t, selecetedObs])

        f_Gibbs_seq_plot[[j]] <- local({
          j <- j
          s <- s
          s2 <- s
          t <- t
          if (p_joint == 0) {
            s2 <- s2 + 1
          }
          if (s %in% 1:p_joint & p_joint != 0) {
            title_y <- paste0("", "p = J", s, ", T =", t)
          } else {
            title_y <- paste0("", "p = I", s2 - 1, ", T =", t)
          }
          ggplot2::ggplot(data = f_Gibbs, ggplot2::aes(x = Gibbs_draws, y = Gibbs_f, color = "Gibbs draws")) +
            ggplot2::geom_line() +
            ggplot2::scale_color_manual(values = c("darkgrey")) +
            ggplot2::ylab(title_y) +
            ggplot2::labs(color = "", title = paste("Mean(Gibbs) =", round(mean(Gibbs$f[s, t, selecetedObs]), 2), ", acf1 =", round(acf(f_Gibbs$Gibbs_f, plot = F)[[1]][2], 2))) +
            ggplot2::theme_bw()
        })
        j <- j + 1
      }
    }

    if (sampleV) {
      # startV <- Gibbs$initials$Vhat[1,1,1]
      V_plot <- vector("list", 2 * N * npara)
      V_data <- Gibbs$V[, selecetedObs]
      BayEstV <- apply(V_data, 1, mean)
      j <- 1
      for (i in 1:(N * npara)) {
        V_Gibbs <- dplyr::tibble(Gibbs_draws = 1:dim(V_data)[2], Gibbs_V = V_data[i, ])
        V_plot[[j]] <- local({
          i <- i
          if (is.null(nameMat)) {
            ylab_text <- paste0("Var", i)
          } else {
            ylab_text <- paste("Var", nameMat[i, 1], nameMat[i, 2])
          }
          # start_val <- startV
          ggplot2::ggplot(data = V_Gibbs, mapping = ggplot2::aes(x = Gibbs_draws, y = Gibbs_V, color = "Gibbs draws")) +
            ggplot2::geom_line() +
            ggplot2::geom_hline(ggplot2::aes(yintercept = BayEstV[i], color = "Bay. Est.")) +
            # ggplot2::labs(color="",title=paste("Est =",round(BayEstV[i],2),", Start =", start_val,", acf1 =", round(acf(V_Gibbs$Gibbs_V,plot=F)[[1]][2],2)) )+
            ggplot2::labs(color = "", title = paste("Est =", round(BayEstV[i], 2), ", acf1 =", round(acf(V_Gibbs$Gibbs_V, plot = F)[[1]][2], 2))) +
            ggplot2::ylab(ylab_text) +
            ggplot2::theme_bw()
        })

        j <- j + 1
        V_plot[[j]] <- local({
          i <- i
          if (is.null(nameMat)) {
            ylab_text <- paste("ACF Var", i)
          } else {
            ylab_text <- paste("ACF", nameMat_ext[i, 1], nameMat_ext[i, 2])
          }
          bacf <- acf(V_Gibbs$Gibbs_V, plot = FALSE)
          bacfdf <- with(bacf, data.frame(lag, acf))
          ggplot2::ggplot(data = bacfdf, mapping = ggplot2::aes(x = lag, y = acf)) +
            ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
            ggplot2::geom_segment(mapping = ggplot2::aes(xend = lag, yend = 0)) +
            ggplot2::ylab(ylab_text)
        })
        j <- j + 1
      }

      return(list(y_plot = y_plot, B_plot = B_plot, D_plot, D_plot, f_plot = f_plot, f_Gibbs_seq_plot = f_Gibbs_seq_plot, Bf_Gibbs_seq_plot = Bf_Gibbs_seq_plot, V_plot = V_plot))
    }


    if (sampleA) {
      startA <- Gibbs$initials$A
      startA_lt <- c(startA[lower.tri(startA, diag = T)])

      if (countryA) {
        A_data <- apply(Gibbs$A[, , selecetedObs, ], 3, function(x) {
          x[lower.tri(Gibbs$A[, , 1, 1], diag = T)]
        })
        A_plot <- vector("list", N * npara * (npara + 1) / 2)
      } else {
        A_data <- apply(Gibbs$A[, , selecetedObs], 3, function(x) {
          x[lower.tri(Gibbs$A[, , 1], diag = T)]
        })
        A_plot <- vector("list", npara * (npara + 1) / 2)
      }

      BayEstA <- apply(A_data, 1, mean)
      name_grid <- expand.grid(1:3, 1:3)
      name_grid <- name_grid[name_grid[, 1] >= name_grid[, 2], ]
      A_names <- as.character(apply(name_grid, 1, \(x) paste0(x[1], x[2])))
      if (countryA) {
        A_names <- rep(A_names, N)
        startA_lt <- rep(startA_lt, N)
      }

      for (i in 1:length(A_plot)) {
        A_Gibbs <- dplyr::tibble(Gibbs_draws = 1:dim(A_data)[2], Gibbs_A = A_data[i, ])
        A_plot[[i]] <- local({
          i <- i
          name <- A_names[i]
          start_val <- startA_lt[i]
          ggplot2::ggplot(data = A_Gibbs, mapping = ggplot2::aes(x = Gibbs_draws, y = Gibbs_A, color = "Gibbs draws")) +
            ggplot2::geom_line() +
            ggplot2::geom_hline(ggplot2::aes(yintercept = BayEstA[i], color = "Bay. Est.")) +
            ggplot2::labs(color = "", title = paste("Est =", round(BayEstA[i], 2), "Start =", start_val, ", acf1 =", round(acf(A_Gibbs$Gibbs_A, plot = F)[[1]][2], 2))) +
            ggplot2::ylab(paste0("A", name)) +
            ggplot2::theme_bw()
        })
      }
      return(list(y_plot = y_plot, B_plot = B_plot, D_plot = D_plot, f_plot = f_plot, f_Gibbs_seq_plot = f_Gibbs_seq_plot, Bf_Gibbs_seq_plot = Bf_Gibbs_seq_plot, A_plot = A_plot, Mu_plot = Mu_plot, Omega_plot = Omega_plot))
    }
  } else {
    B_plot <- NULL
    D_plot <- NULL
    f_plot <- NULL
    f_Gibbs_seq_plot <- NULL
    Bf_Gibbs_seq_plot <- NULL
    Mu_plot <- Mu_plot <- NULL
    Omega_plot <- NULL
  }

  return(list(y_plot = y_plot, B_plot = B_plot, D_plot = D_plot, f_plot = f_plot, f_Gibbs_seq_plot = f_Gibbs_seq_plot, Bf_Gibbs_seq_plot = Bf_Gibbs_seq_plot, Mu_plot = Mu_plot, Omega_plot = Omega_plot))
}
save_zoom <- function(plot_list, ncol, as.table = T, top = "", file = "plot", path, width = 1536 / 72, height = 814 / 72, dpi = 72) {
  plotArr <- do.call(eval(parse(text = "gridExtra::arrangeGrob")), c(plot_list, ncol = ncol, as.table = as.table, top = top))
  if (!missing(path)) {
    ggplot2::ggsave(paste0(path, paste0("/", file, ".pdf")), plotArr, width = width, height = height, dpi = dpi)
  } else {
    ggplot2::ggsave(paste0(file), plotArr, width = width, height = height, dpi = dpi)
  }
}
