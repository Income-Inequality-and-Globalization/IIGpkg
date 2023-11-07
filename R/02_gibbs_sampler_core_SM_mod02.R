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
#'    `VhatDiagScale = TRUE` ist muss auch `VdiagEst = TRUE` sein.
#' @param VhatDiagScale_start Startwert fuer die Skalierungsfaktoren
#' @param VdiagEst (wahrscheinlich irrelevant) = TRUE, wenn die diagonale Messfehlervarianz geschaetzt werden soll. (Fuer uns irrelevant, wir wollen die GMM-Varianz nutzen)
#' @param alpha0 Prior-Shape-Paramter fuer Inverse-Gamma(!), wenn VdiagEst = TRUE
#' @param beta0 Prior-Scale-Paramter fuer Inverse-Gamma(!), wenn VdiagEst = TRUE
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
                       B0,
                       Omega0,
                       selectR,
                       selectC,
                       Vhat,
                       incObsOld = 100000,
                       incObsNew = 100000,
                       covScale,
                       VhatDiagScale,
                       VhatDiagScale_start = NULL,
                       VdiagEst,
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
  initials <- as.list(environment()) # speichert die Initialisierung
  # DEBUG_ITER <- 464
  Vhat             <- set_scale_Vhat(Vhat, incObsOld, incObsNew)
  store_count      <- 0 # zaehlt wie oft bereits zwischengespeichert wurde
  SIMULATE_FACTORS <- set_simulate_factors(fPost)
  w_reg_spec       <- set_reg_specification(wRegSpec)
  msg_error_kf     <- NULL

  # browser()
  num_y   <- npara
  N_num_y <- dim(B)[1] # N * number of components in y
  num_fac <- dim(B)[2] # number of factors
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
  invOmega0 <- solve(Omega0)
  invOmega0_B0_D0 <- compute_invOmega0_B0_D0(invOmega0, B0, D0) 

  # Container fuer Gibbs-Zuege fuer f, B, D, A, V
  fSTORE <- set_f_out(num_fac, TT, itermax)
  BSTORE <- set_B_out(N_num_y, num_fac, itermax)
  DSTORE <- set_D_out(N_num_y, NN, nreg, itermax)
  if (is.null(DSTORE)) {
    D <- matrix(0, nrow = N_num_y)
    wReg <- matrix(0, ncol = TT)
  }
  if (nreg == 0) D0 <- NULL
  ASTORE <- set_A_out(countryA, num_y, itermax, NN) 
  VSTORE <- set_V_out(N_num_y, itermax)
  V_tmp  <- set_V_tmp(VdiagEst, VhatDiagScale,
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
                                              num_reg = nreg,
                                              i = 2)
  upper <- tmp_ir$upper
  lower <- tmp_ir$lower

  # storePath Anpassung
  store_paths <- set_store_path_subdir(
    storePath, VdiagEst, sampleA,
    initials, Vstart, covScale, njointfac,
    w_reg_spec, Omega0, incObsNew
  )
  storePath_adj <- store_paths$store_path_adj
  storePath_rds <- store_paths$store_path_rds
  storePath_omg <- store_paths$store_path_omg
  storePath_kfe <- store_paths$store_path_kfe
  ##############################################################################
  ####################### GIBBS sampler Iteration START ########################
  ##############################################################################
  # Wird nicht mehr benoetigt, nicht auskommentiert, da es unten auch noch drin
  # steht. War vom alten Sampler, um die identifizierenden Restriktionen
  # einzuhalten:
  block_count <- 0 
  # Gibbs-Iteration:
  for (iter in seq_len(itermax)) {
    # if (iter == DEBUG_ITER) browser()
    # Erweiterung das Vhat-Arrays um die Adjustmentmatrix A
    VhatArray_A <- compute_V_hat_array_A(VhatSqrt, A,
                                         countryA, A_countryArray,
                                         NN, TT, NN_TT, num_y)
    #### GIBBS PART: Sampling of the latent factors (FFBS)
    if (SIMULATE_FACTORS) {
      fPost <- compute_FFBS(yObs = yObs, uReg = uReg, wReg = wReg,
                            num_fac = num_fac, num_y = num_y, N_num_y = N_num_y,
                            TT = TT, NN = NN,
                            VhatArray_A = VhatArray_A, Phi = Phi,
                            B = NULL, C = B, D = D, Q = Q,
                            initX = initX, initU = initU, initP = initP,
                            PDSTORE = FALSE,
                            try_catch_errors = NULL) #,
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
      fSTORE[, , iter] <- fPost
    } else {
      fSTORE[, , iter] <- fPost
    }
    ##### GIBBS PART: Sampling of loadings (on latent factors) and partial effects (on regressors)
    # Ergebnis-Matrix fuer die gemeinsamen Ladungen (wird spaeter befuellt)
    Bjoint <- matrix(rep(0, N_num_y * njointfac), ncol = njointfac)
    # Ergebnis-Matrix fuer die idiosynkratischen Ladungen (wird spaeter befuellt)
    Bidio <- matrix(rep(0, N_num_y), ncol = NN)

    for (i in 1:NN) {
      # Vhat Array for cross-sectional unit (country) i
      Viarray <- VhatArray_A[, , (1 + (i - 1) * TT):(i * TT)]
      # yObs for cross-sectional unit (country) i
      yiObs <- yObs[(1 + num_y * (i - 1)):(num_y * i), ]
      # availableObs <- which(!is.na(yObs[1 + num_y * (i-1),]))
      availableObs <- which(availableObs_crossSection[i, ])

      Omega1 <- compute_Omega1(invOmega0 = invOmega0,
                               availableObs = availableObs,
                               selectR = selectR,
                               num_y = num_y,
                               nreg = nreg,
                               njointfac = njointfac, 
                               i = i,
                               fPost = fPost,
                               wReg = wReg,
                               Viarray = Viarray,
                               type = type,
                               storePath_adj = storePath_adj,
                               storePath_omg = storePath_omg,
                               store_count = store_count,
                               initials = initials,
                               incObsNew = incObsNew,
                               iter = iter)
      beta1 <- compute_B_mean(Omega = Omega1,
                              invOmega0_B0_D0[[i]],
                              availableObs = availableObs,
                              selectR = selectR,
                              num_y = num_y,
                              nreg = nreg, 
                              njointfac = njointfac,
                              i = i,
                              fPost = fPost,
                              wReg = wReg,
                              yiObs = yiObs, 
                              Viarray = Viarray,
                              type = type)
      Sigma <- compute_Sigma_adjust(Omega1)
      # Sampling der Ladungen bzw. partiellen Effekte
      B_D_samp <- sample_B_D(mean_B_full = beta1, Sigma,
                             upper = upper, lower = lower,
                             num_jnt_fac = njointfac, num_y = num_y)
      Dsamp <- B_D_samp$Dsamp
      Bsamp_jnt <- B_D_samp$Bsamp_jnt
      Bsamp_idi <- B_D_samp$Bsamp_idi

      if (nreg != 0) {
        DSTORE[(1 + num_y * (i - 1)):(i * num_y),
               (1 + nreg * (i - 1)):(i * nreg), iter] <- Dsamp
      }
      if (njointfac != 0) {
        Bjoint[((i - 1) * num_y + 1):(i * num_y), ] <- Bsamp_jnt
      }
      Bidio[, i] <- Bsamp_idi
    }
    B <- get_B(Bjoint, Bidio, num_fac_jnt = njointfac,
               NN = NN, num_y = num_y, type = type)
    BSTORE[, , iter] <- B
    if (nreg != 0) D <- DSTORE[, , iter]
    # if (iter == DEBUG_ITER) browser()
    ############################################################################
    ######### GIBBS PART: Sampling VCOV matrix or Adjustment-Matrix A ##########
    ############################################################################
    u <- compute_residuals(yObs, fPost, B, D, wReg) # uSTORE[,, iter] <- u
    if (VdiagEst) {
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
        A <- sample_A(countryA, diagA, scaleA,
                      u, VhatSqrt,
                      NN, TT, num_y, TT_num_y,
                      NN_TT_avail, availableObs_crossSection,
                      shape0, rate0, Psi0, nu0)
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
      store_mcmc(VdiagEst, initials,
                 fSTORE, BSTORE, DSTORE, VSTORE, ASTORE,
                 block_count, msg_error_kf, storePath_rds)
        
    }
    print(iter)
  }
  return(list(f = fSTORE, B = BSTORE, D = DSTORE, A = ASTORE, V = VSTORE, 
              blockCount = block_count, errorMsg = msg_error_kf,
              initials = initials))
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
