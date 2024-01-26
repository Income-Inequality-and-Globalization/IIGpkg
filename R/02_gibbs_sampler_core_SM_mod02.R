#' Title
#'
#' @param itermax Number of Gibbs iterations
#' @param identmax Nicht mehr noetig (war fuer altes Sampling proceduere - so lange samplen bis Parameterrestriktionen eingehalten)
#' @param npara Anzahl der Parameter (bei uns immer 3)
#' @param nreg Anzahl der  Regressoren (1 bis 3)
#' @param njointfac Anzahl der gemeinsamen Faktoren (zuletzt immer 0)
#' @param yObs Daten (a, q, mu)
#' @param c_ Konstante der Transition-equation (keine Konstante im Modell - daher nicht relevant)
#' @param D Startwert fuer die Koeffizienten der Regressorn (stacked: npara*N x nreg*NN), N = Anzahl Laender
#' @param D0 Prior Erwwartungswert fuer die Koeffizienten der Makro- Regressoren (bei mir immer gleich dem Starwert, da ich beides mit makeDstart erzeuge. Aber nicht stacked, sondern npara x nreg x N Array)
#' @param B Startwerte fuer die Ladungen auf die latenten Faktoren (stacked: npara *N x #Faktoren)
#' @param Phi Koeffizentenmatrix des State-Prozesses (bei unser immer diag(#Faktoren), random walks)
#' @param Q Felhervarianz des State-Prozesses (bi uns immer auf diag(#Faktoren) gesetzt)
#' @param type Modelltyp: Auswahls aus "allidio","countryidio", "countryidio_nomu". Bei uns immer "allidio", d.h. jeder Paramter (jeden Landes) bekommt eigenen idiosynkr. Faktor. Hatte auch mal ausprobiert, dass z.B. nur jedes Land einen bekommt.
#' @param initX Initialsierung der States im KF. Bei uns immer rep(0, #Faktoren)
#' @param initP Intitalisierung der State-Varianz im KF. Bei uns immer diag(0, #Faktoren)
#' @param initU Initialisierung der Regressoren in der State-Prozess. Für uns irrelevant, darf im KF aber nciht fehlen und wird daher = 0.
#' @param wRegSpec Ein Skalar, der angibt, welche Regressoren hinzugefuegt werden. Bei drei Regressoren gibt es 7 Möglichkeiten, also ist wRegSpec eine Zahl von 1-7. Die Moeglichkeiten werden in Gibbs2_SM_SA_cluster berechnet. Es gilt: 1 = "cpi_change", 2 = "unemployment", 3 = "gdp_ppp", 4 = "cpi_change" + "unemployment", 5 = "cpi_change" + "gdp_ppp", 6 = "unemployment" + "gdp_ppp", 7 = "cpi_change" +"unemployment" +"gdp_ppp"
#' @param wReg Regressoren für die Messgleichung, also die Makros: nreg * N x T (muss natuerlich zu wRegSpec passen)
#' @param uReg Siehe initU. Fuer uns irrelevant, wird gleich NULL gesetzt.
#' @param B0 Prior Erwwartungswert fuer die Ladungen auf die latenten Fakoren (bei mir immer gleich dem Starwert, da ich beides mit makeBstart erzeuge. Aber nicht stacked, sondern npara x (npara + njointfac) x N Array bei type = "allidio")
#' @param Omega0 Prior Varianz fuer die Ladungen auf die latenten Fakoren und die Makro-Regressor-Koeffizienten.
#' @param selectR Selektionsmatrix: Waehlt Elemente aus vec(B_i), die gesamplet werden sollen. Auf den anderen Elementen liegen Nullrestriktionen
#' @param selectC Es wird selectR * vec(B_i) + selectC gesamplet. selectC ist bei uns immer ein Nullvektor, da wir keine affin-lineare Transformation benoetigen.
#' @param Vhat QML(bzw. GMM) geschaetzte Fehlervarianz der Vorschaetzung. Ein npara x npara x N * T Array. Die ersten Elemente [1:npara, 1:npara, 1:T] sind die Varianzen des ersten Landes ueber die Zeit, [1:npara, 1:npara, (T+1):(2T)] sind die Varianzen des zweiten Landes ueber die Zeit, usw.
#' @param incObsOld Anzahl der Beobachtungen (der Rohdaten, also individuelle Einkommen), die wir in der GMM-Schaetzung angenommen haben. Bei mir 100000. Das ist zu viel. Daher incObsNew.
#' @param incObsNew Gewuenschte / Korrigierte Anzahl der Beochbachtungen in der GMM-Schaetzung. Setze bspw. gleich 5000 oder 10000. Das ist realistischer als 100000.
#' @param covScale Damit kann man die Kovarianzen aus Vhat skalieren. covScale = 0.5 heißt, dass die die Korrelationen(!) halbiert werden. War mal ein Experiment. covscale muss letztlich = 1 sein.
#' @param VhatDiagScale (wahrscheinlich irrelevant) `= TRUE`, wenn fuer jedes
#'    Diagonalelement der GMM- varianz ein Skalierungsfaktor geschaetzt werden 
#'    soll. `Diag_VCOV = TRUE` muss in der top-level Function gesetzt werden,
#'    damit die Nebendiagonalelemente der GMM-Varianz = 0 gesetzt werden. Wenn
#'    `VhatDiagScale = TRUE` ist muss auch `sampleV = TRUE` sein.
#' @param VhatDiagScale_start Startwert fuer die Skalierungsfaktoren
#' @param sampleV (wahrscheinlich irrelevant) = TRUE, wenn die diagonale Messfehlervarianz geschaetzt werden soll. (Fuer uns irrelevant, wir wollen die GMM-Varianz nutzen)
#' @param alpha0 Prior-Shape-Paramter fuer Inverse-Gamma(!), wenn sampleV = TRUE
#' @param beta0 Prior-Scale-Paramter fuer Inverse-Gamma(!), wenn sampleV = TRUE
#' @param countryA = TRUE, wenn fuer jedes Land eine eigene Adjustmentmatrix (der GMM-Varianz) A geschaetzt werden soll
#' @param A Startwert fuer A, bei mir immer diag(npara)
#' @param scaleA =TRUE, wenn die Matrix A nur ein Skalar sein soll.
#' @param diagA = TRUE, wenn die Matrix A eine Diagonalmatrix sein soll
#' @param Psi0 Prior-Paramter der Inverse-Wishart (wenn A eine vollbesetzte Matrix ist)
#' @param nu0 Prior-Paramter der Inverse-Wishart (wenn A eine vollbesetzte Matrix ist)
#' @param shape0 Prior-Shape-Paramter fuer Inverse-Gamma(!), wenn diagA = TRUE oder scaleA = TRUE
#' @param rate0 Prior-Scale-Paramter fuer Inverse-Gamma(!), wenn diagA = TRUE oder scaleA = TRUE
#' @param storePath Speicherpfad des .rds Outputs. Muss nur der Oberordner sein. Es werden automatisch Unterordner mit entpsrechenden Ordnernamen angelegt. Setze = "none", wenn nichts gespeichert werden soll.
#' @param fPost Kann angegeben werden, wenn FFBS ausgelassen werden soll und man gegeben fester Faktoren samplen möchte. fPost muss dann #Faktoren x T sein. Wenn es fehlt, wird FFBS durchgefuehrt.
#' @param sampleA = TRUE, wenn A gesamplet werden soll (bei uns immer der Fall). Bei sampleA = FALSE, muss trotzdem das A Argument = diag(npara) gesetzt werden.
#' @param identification (Kann ignoriert werden) Boolscher Wert der dazu diente, um identifizierende Restriktionen zu ignorieren. Jetzt werde diese immer eingehalten.
#' @param storeUnit Zahl. Alle storeUnit-Iterationen werden die Zwischenergebnisse gespeichert.
#'
#' @return core Gibbs sampler output
#' @export
GibbsSSM_2 <- function(itermax = 15000,
                       identmax = 5000,
                       npara,
                       nreg,
                       njointfac,
                       yObs,
                       c_,
                       D,
                       D0,
                       B,
                       Phi,
                       Q,
                       type = "allidio",
                       initX,
                       initP,
                       initU,
                       wRegSpec,
                       wReg,
                       uReg = NULL,
                       prior_list,
                       selectR,
                       selectC,
                       Vhat,
                       incObsOld = 100000,
                       incObsNew = 100000,
                       covScale,
                       VhatDiagScale,
                       VhatDiagScale_start = NULL,
                       sampleV,
                       alpha0,
                       beta0,
                       countryA,
                       A,
                       scaleA,
                       diagA,
                       Psi0,
                       nu0,
                       shape0,
                       rate0,
                       storePath = "none",
                       fPost,
                       sampleA,
                       identification = TRUE,
                       storeUnit = 10000) {
  # Vhat muss (npara, npara, N * T) Array sein. Start: Alle Zeitpunkte fuer das erste Land.
  # Funktioniert noch nicht fuer npara < njointfac (Geweke)
  # Modell ohne Konstante
  # initials <- get_initials(as.list(environment())) # speichert die Initialisierung
  B0       <- prior_list$priors$B0
  Omega0   <- prior_list$priors$Omega0
  sampleH  <- FALSE
  if (!is.null(prior_list$hyperpriors)) {
    mu_b0    <- prior_list$hyperpriors$mu_b0
    Sigma_b0 <- prior_list$hyperpriors$Sigma_b0
    alpha_b0 <- prior_list$hyperpriors$alpha_b0
    beta_b0  <- prior_list$hyperpriors$beta_b0
    sampleH  <- TRUE
  }
  initials <- as.list(environment())
  DEBUG_ITER <- 103
  Vhat             <- set_scale_Vhat(Vhat, incObsOld, incObsNew)
  store_count      <- 0 # zaehlt wie oft bereits zwischengespeichert wurde
  SIMULATE_FACTORS <- set_simulate_factors(fPost)
  w_reg_spec       <- set_reg_specification(wRegSpec)
  msg_error_kf     <- NULL

  num_y   <- npara
  N_num_y <- dim(B)[1] # N * number of components in y
  num_fac <- dim(B)[2] # number of factors
  num_fac_jnt <- njointfac
  num_fac_idi <- num_fac - num_fac_jnt
  num_par_all <- N_num_y + num_fac_jnt
  # Anzahl der Laender:
  NN       <- N_num_y / num_y
  # Anzahl der Zeitpunkte:
  TT       <- dim(yObs)[2]
  # Laender mal Zeitpunkte
  NN_TT    <- NN * TT
  # Zeitpunkte mal y-measurement components:
  TT_num_y <- TT * num_y 
  # N *TT ohne die fehlenden Beobachtungen:
  NN_TT_avail <- sum(!is.na(yObs)) / num_y
  # Anzahl der Beobachtungen (ohne die NAs) fuer jedes Land:
  availableObs_crossSection <- t(
    apply(yObs[seq(1, N_num_y, 3), ],  1, \(x) !is.na(x))
  ) 

  Omega0_legacy <- Omega0
  if (sampleH) {
    invOmega0 <- solve(Omega0)
    # invOmega0_B0_D0 <- NULL
  } else {
    invOmega0 <- solve(Omega0)
    # invOmega0_B0_D0 <- compute_invOmega0_B0_D0(invOmega0, B0, D0)
    # Omega0 <- NULL
    # B0 <- NULL
    # D0 <- NULL
  }

  # Container fuer Gibbs-Zuege fuer f, B, D, A, V
  fSTORE  <- set_f_out(num_fac, TT, itermax)
  BSTORE  <- set_B_out(N_num_y, num_fac, itermax)
  DSTORES <- set_D_out(N_num_y, NN, nreg, itermax)
  DSTORE  <- DSTORES$DSTORE
  DiSTORE <- DSTORES$DiSTORE
  # uSTORE <- array(0, dim = c(N_num_y, TT, itermax) )
  uSTORE <- NULL
  BD0STORE    <- matrix(0, nrow = sum(selectR), ncol = itermax)
  Omega0STORE <- matrix(0, nrow = sum(selectR), ncol = itermax)

  # Availability of regressors
  if (is.null(DSTORE)) {W_REG_AVAIL <- TRUE} else {W_REG_AVAIL <- TRUE}
  if (W_REG_AVAIL) {
    id_w_reg <- get_id_wreg(nreg, NN)
    w_reg_info = list(w_reg = wReg, 
                      id_reg = id_w_reg,
                      nregs = nreg)
  } else {
    D <- matrix(0, nrow = N_num_y)
    wReg <- matrix(0, ncol = TT)
    D0 <- NULL
    w_reg_info <- NULL
  }

  ASTORE <- set_A_out(countryA, num_y, itermax, NN) 
  VSTORE <- set_V_out(sampleV, N_num_y, itermax)
  V_tmp  <- set_V_tmp(sampleV, VhatDiagScale,
                      Vhat, VhatDiagScale_start,
                      NN, TT, NN_TT, N_num_y, num_y)
  VhatArrayBdiagByTimeFix <- V_tmp$VhatArrayBdiagByTimeFix
  VhatFix                 <- V_tmp$VhatFix
  Vstart                  <- V_tmp$Vstart
  Vhat                    <- V_tmp$Vhat

  VhatSqrt <- compute_mat_square_root_V(Vhat, num_y, NN_TT)
  A_countryArray <- set_A_country_array(countryA, A, NN, num_y)
  # pb <- winProgressBar(title="Progress", label="0% done", min=0, max=100, initial=0)
  # ident_block <- FALSE # wird nicht mehr benoetigt, nicht auskommentiert, da
  # es unten auch noch drin steht. War vom alten Sampler, um die
  # identifizierenden Restriktionen einzuhalten.
  tmp_ir <-  get_identificiation_restrictions(type = type,
                                              num_joint_fac = njointfac,
                                              num_y = num_y,
                                              num_reg = nreg)
  upper <- tmp_ir$upper
  lower <- tmp_ir$lower

  # storePath Anpassung
  store_paths <- set_store_path_subdir(
    storePath, sampleV, sampleA, sampleH,
    initials, Vstart, covScale, njointfac,
    w_reg_spec, Omega0_legacy, incObsNew
  )
  storePath_adj <- store_paths$store_path_adj
  storePath_rds <- store_paths$store_path_rds
  storePath_omg <- store_paths$store_path_omg
  storePath_kfe <- store_paths$store_path_kfe

  id_f     <- get_id_fpost(num_fac_jnt, num_y, NN, type)
  ##############################################################################
  ####################### GIBBS sampler Iteration START ########################
  ##############################################################################
  # Wird nicht mehr benoetigt, nicht auskommentiert, da es unten auch noch drin
  # steht. War vom alten Sampler, um die identifizierenden Restriktionen
  # einzuhalten:
  block_count <- 0 
  # Gibbs-Iteration:
  # load("kf_test2.RData")
  # load("kf_test.RData")
  # browser()
  SIMULATE_FACTORS <- TRUE
  # fPost_true <- all_true_vals$f_simul
  fPost_true <- NULL
  # A <- all_true_vals$A_means
  # A <- diag(nrow(A))
  # B <- all_true_vals$B_means
  B[, 1] <- 1; diag(B[, 2:31]) <- 1; #B[which(B != 0)] <- 0.5 # 
  # D <- all_true_vals$D_means
  # D[which(D != 0)] <- 0.5
  for (iter in seq_len(itermax)) {
    # if (iter == DEBUG_ITER) browser()
    # Erweiterung das Vhat-Arrays um die Adjustmentmatrix A
    VhatArray_A <- compute_V_hat_array_A(VhatSqrt, A,
                                         countryA, A_countryArray,
                                         NN, TT, NN_TT, num_y)
    #### GIBBS PART: Sampling of the latent factors (FFBS)
    if (SIMULATE_FACTORS) {
      # R <- bdiagByTime(
      #   VhatArray_A = VhatArray_A,
      #   npara = num_y,
      #   N = NN,
      #   TT = TT,
      #   Nnpara = N_num_y
      # )
      # fPost <- compute_FFBS2(yObs = yObs, uReg = uReg, wReg = wReg,
      #                        num_fac = num_fac, num_y = num_y, N_num_y = N_num_y,
      #                        TT = TT, NN = NN,
      #                        R = R, Phi = Phi,
      #                        B = NULL, C = B, D = D, Q = Q,
      #                        initX = initX, initU = initU, initP = initP,
      #                        PDSTORE = FALSE,
      #                        try_catch_errors = NULL) #,
      # try_catch_errors =
      #   list(storePath_adj = storePath_adj,
      #        store_count = store_count,
      #        fSTORE = fSTORE,
      #        BSTORE = BSTORE,
      #        ASTORE = ASTORE,
      #        storePath_kfe = storePath_kfe,
      #        block_count = block_count,
      #        initials = initials,
      #        iter = iter))
      fPost <- ffbs(yObs, wReg,
                    num_fac, N_num_y, TT,
                    initX, initP,
                    Phi, B, D, Q, VhatArray_A,
                    FALSE)
      # fPost[2:31, ] <- fPost_true[2:31, ]
      # fPost[1, ] <- fPost_true[1, ]
      fSTORE[, , iter] <- fPost
      plot_mcmc_ffbs(f_stored = fSTORE,
                     mm = iter,
                     f_trues = fPost_true,
                     options_f = list(
                             grid = c(4, 4),
                             id_f = 1:16),
                     options_p = list(
                             break_mm = 10000,
                             pause = 0.1))
    } else {
      fPost <- fPost_true
      fSTORE[, , iter] <- fPost_true
    }
    ##### GIBBS PART: Sampling of loadings (on latent factors) and partial
    ##### effects (on regressors)
    B_post_all <- sample_B_full_cpp(yObs, availableObs_crossSection,
                                    fPost, VhatArray_A, w_reg_info,
                                    Omega0 = Omega0, B0 = B0, D0 = D0,
                                    id_f, selectR, lower, upper,
                                    NN, TT, N_num_y, num_y, num_fac_jnt, num_par_all,
                                    type, 
                                    B_trues = all_true_vals$B_means,
                                    D_trues = all_true_vals$D_means)
                                # if (sampleH),)
    # B_post_all <- sample_B_full(yObs, availableObs_crossSection,
    #                             fPost, VhatArray_A, w_reg_info,
    #                             Omega0 = Omega0, B0 = B0, D0 = D0,
    #                             id_f, selectR, lower, upper,
    #                             NN, TT, N_num_y, num_y, num_fac_jnt, num_par_all,
    #                             type)
    DSTORE[, , iter] <- B_post_all$Dregs
    DiSTORE          <- B_post_all$Dregs_i # if (sampleH)
    BSTORE[, , iter] <- B_post_all$Bfacs
    B   <- B_post_all$Bfacs
    B_i <- B_post_all$Bfacs_i 
    if (nreg != 0) D <- DSTORE[, , iter]
    #### GIBBS Part : Sampling of B0, D0, and Omega0
    if (sampleH) {
      BD_0 <- mu_sampler(NN = NN, npara = num_y, nreg = nreg, 
                         njointfac = num_fac_jnt,
                         B_i = B_i, 
                         D_i = DiSTORE,
                         invOmega0 = invOmega0, mu_b0 = mu_b0, Sigma_b0 = Sigma_b0, 
                         selectR = selectR, type = type)
      B0 <- BD_0$B0
      D0 <- BD_0$D0
      BD0STORE[, iter] <- selectR %*% c(B0[,,1], D0[,,1])

      Omega0 <- omega_sampler(NN = NN, 
                              B_i = B_i, D_i = DiSTORE,
                              B0 = B0, D0 = D0, 
                              alpha_b0 = alpha_b0, beta_b0 = beta_b0, 
                              selectR = selectR)
      Omega0STORE[, iter] <- diag(Omega0)
      invOmega0 <- solve(Omega0)
    }
    # if (iter == DEBUG_ITER) browser()
    ############################################################################
    ######### GIBBS PART: Sampling VCOV matrix or Adjustment-Matrix A ##########
    ############################################################################
    u <- compute_residuals(yObs, fPost, B, D, wReg) # uSTORE[,, iter] <- u
    if (sampleV) {
      tmp_V <- sample_V(VhatDiagScale, VhatArrayBdiagByTimeFix, VhatFix, u,
                        TT, NN_TT, num_y, availableObs_crossSection,
                        alpha0, beta0)
      VhatSqrt       <- tmp_V$V_hat_sqrt
      VSTORE[, iter] <- tmp_V$V_diag
    }
    if (sampleA) {
      if (countryA) {
        A_countryArray <- sample_A(countryA, diagA, scaleA,
                                   u, VhatSqrt,
                                   NN, TT, num_y, TT_num_y,
                                   NN_TT_avail, availableObs_crossSection,
                                   shape0, rate0, Psi0, nu0)
        ASTORE[, , iter, ] <- A_countryArray
      } else {
        # browser()
        A <- sample_A(countryA, diagA, scaleA,
                      u, VhatSqrt,
                      NN, TT, num_y, TT_num_y,
                      NN_TT_avail, availableObs_crossSection,
                      shape0, rate0, Psi0, nu0)
        # A <- all_true_vals$A_means
        ASTORE[, , iter] <- A
      }
    }
    CHECK_STORE <- get_check_store(store_path = storePath,
                                   mcmc_iter = iter,
                                   mcmc_itermax = itermax,
                                   store_unit = storeUnit)
    if (CHECK_STORE) {
      store_count <- store_count + 1
      if (store_count == 1) dir.create(storePath_adj, recursive = TRUE)
      out_list <- get_outlist_Gibbs_sampler(fSTORE, BSTORE, DSTORE,
                                            ASTORE, VSTORE, uSTORE,
                                            BD0STORE, Omega0STORE,
                                            block_count,
                                            msg_error_kf,
                                            initials)
      store_mcmc(out_list, storePath_rds)
    }
    print(iter)
  }
  out_list <- get_outlist_Gibbs_sampler(fSTORE, BSTORE, DSTORE,
                                        ASTORE, VSTORE, uSTORE,
                                        BD0STORE, Omega0STORE,
                                        block_count,
                                        msg_error_kf,
                                        initials)
  out_list <- new_GibbsOutputIIG(out_list)
  return(new_GibbsOutputIIG(out_list))
}
############################################################################
## GIBBS sampler Iteration ENDE
############################################################################
# info <- sprintf("%d%% done", round((iter/itermax)*100))
# setWinProgressBar(pb, iter/itermax*100, label=info)
############################################################################
## GIBBS sampler SPEICHERN
############################################################################
# if(iter %% 10000 == 0){print(iter)}
# close(pb)
# RETURN WITH uSTOREÖ
# return(list(f = fSTORE, B = BSTORE, D = DSTORE, A = ASTORE, V = VSTORE,
# u = uSTORE, blockCount = blockCount, errorMsg = errorMsg, initials = initials))
# valid <- FALSE
# ident_control <- 1
# while(!valid){
#   if(ident_control == identmax + 1){iter <- iter - 1; ident_block <- T; blockCount <- blockCount +1 ; break}
