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

  # Container fuer Gibbs-Zuege fuer f, B, D, A, V
  fSTORE <- set_f_out(num_fac, TT, itermax)
  BSTORE <- set_B_out(N_num_y, num_fac, itermax)
  DSTORE <- set_D_out(N_num_y, NN, nreg, itermax)
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
  # storePath Anpassung
  storePath_adj <- set_store_path_subdir(
    storePath, VdiagEst, sampleA,
    initials, covScale, njointfac,
    w_reg_spec, Omega0, incObsNew
  )
  ##############################################################################
  ####################### GIBBS sampler Iteration START ########################
  ##############################################################################
  # Wird nicht mehr benoetigt, nicht auskommentiert, da es unten auch noch drin
  # steht. War vom alten Sampler, um die identifizierenden Restriktionen
  # einzuhalten:
  block_count <- 0 
  # Gibbs-Iteration:
  iter <- 1
  for (iter in seq_len(itermax)) {
    # if (iter == 5) browser()
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
                            storePath_adj = storePath_adj, store_count = store_count)
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


      ## Posterior Momente fuer die Ladungen und die partiellen Effekte
      invOmega0 <- solve(Omega0)
      invOmega1_part2 <- sumffkronV(availableObs, npara = num_y, nreg = nreg, njointfac = njointfac, i = i, fPost = fPost, wReg = wReg, Viarray = Viarray, type = type)
      invOmega1 <- invOmega0 + invOmega1_part2

      beta1_mid <- sumfyV(availableObs, npara = num_y, nreg = nreg, njointfac = njointfac, i = i, fPost = fPost, wReg = wReg, yiObs = yiObs, Viarray = Viarray, type = type)

      # Omega1 <-  tryCatch({solve(invOmega1)},
      #          error = function(e){
      #            #print(invOmega1)
      #            if(storePath != "none"){saveRDS(invOmega1, file =paste0(storePath,"invOmega1_","pj",njointfac,"_B",initials$B0[1,1,1],"_Omega",initials$Omega0[1,1],"_A",initials$A[1,1],"_Psi",initials$Psi0[1,1],"_nu0",initials$nu0,"_",iter))}
      #            #browser()
      #            solve(invOmega1, tol = 0)
      #            })
      #

      Omega1 <- tryCatch(
        {
          solve(selectR %*% invOmega1 %*% t(selectR))
        },
        error = function(e) {
          # print(invOmega1)
          if (storePath != "none") {
            if (store_count == 0) {
              dir.create(storePath_adj, recursive = T)
              store_count <- store_count + 1
            }
            saveRDS(invOmega1, file = paste0(storePath_adj, "/", "invOmega1_", "pj", njointfac, "_B", round(initials$B0[1, 1, 1], 2), "_Omega", initials$Omega0[1, 1], "_A", initials$A[1, 1], "_Psi", initials$Psi0[1, 1], "_nu0", initials$nu0, "_IO", incObsNew, "_", iter))
          }
          # browser()
          solve(selectR %*% invOmega1 %*% t(selectR), tol = 0)
        }
      )



      # beta1 <- Omega1 %*% (c(beta1_mid) + invOmega0 %*%  c(B0[,,i]) )
      beta1 <- Omega1 %*% (selectR %*% (c(beta1_mid) + invOmega0 %*% c(B0[, , i], D0[, , i])))


      # valid <- FALSE
      # ident_control <- 1
      # while(!valid){
      #   if(ident_control == identmax + 1){iter <- iter - 1; ident_block <- T; blockCount <- blockCount +1 ; break}

      ## Identifikationsrestriktionen fuer die Ladungen
      upper <- rep(Inf, 6)
      # Dselect <- diag(6)
      if (type == "allidio") {
        if (njointfac == 0) {
          lower <- c(rep(0, num_y), rep(-Inf, num_y * nreg)) # muss fuer gemeinsamen faktor angepasst werden
          upper <- rep(Inf, num_y + num_y * nreg) # muss fuer gemeinsamen faktor angepasst werden
        } else {
          if (i == 1) {
            lower <- c(0, rep(-Inf, 2), rep(0, 3))
            #   Dselect <- Dselect[-c(2,3),]
          } else {
            lower <- c(rep(-Inf, 3), rep(0, 3))
            #   Dselect <- Dselect[-c(1,2,3),]
          }
        }
      } else if (type == "countryidio") {
        if (i == 1) {
          lower <- c(0, rep(-Inf, 2), 0, rep(-Inf, 2))
          #   Dselect <- Dselect[-c(2,3),]
        } else {
          lower <- c(rep(-Inf, 3), 0, rep(-Inf, 2))
          #   Dselect <- Dselect[-c(1,2,3),]
        }
      } else if (type == "countryidio_nomu") {
        if (i == 1) {
          lower <- c(0, rep(-Inf, 2), 0, -Inf)
          #   Dselect <- Dselect[-c(2,3),]
        } else {
          lower <- c(rep(-Inf, 3), 0, -Inf)
          #   Dselect <- Dselect[-c(1,2,3),]
        }
        upper <- rep(Inf, 5)
      } else {
        stop("No valid model type selected.")
      }
      # muss verallgemeinert werden fuer num_y < njointfac (fuer uns nicht noetig)


      # Bvec <- MASS::mvrnorm(n = 1, mu = selectR %*% beta1 + selectC, Sigma = selectR %*% Omega1 %*% t(selectR))
      # Sigma <- 0.5 * (selectR %*% Omega1 %*% t(selectR)) + 0.5 * t(selectR %*% Omega1 %*% t(selectR))
      Sigma <- 0.5 * Omega1 + 0.5 * t(Omega1)
      # Bvec <- as.numeric(tmvtnorm::rtmvnorm(n = 1, mean = as.numeric(selectR %*% beta1 + selectC), sigma = Sigma,
      #                          lower = lower, upper = upper))

      # Bvec <- tmvnsim::tmvnsim(1,length(upper),lower = lower, upper = upper, means = as.numeric(selectR %*% beta1 + selectC), sigma = Sigma )$samp
      # Bvec <- tmvnsim::tmvnsim(1,length(upper),lower = lower, upper = upper, means = as.numeric(beta1 + selectC), sigma = Sigma )$samp

      # Sampling der Ladungen bzw. partiellen Effekte
      BDsamp <- tmvnsim::tmvnsim(1, length(upper), lower = lower, upper = upper, means = as.numeric(beta1 + selectC), sigma = Sigma)$samp

      Bvec <- BDsamp[1:((njointfac + 1) * num_y)]
      Dvec <- BDsamp[-(1:((njointfac + 1) * num_y))]
      D_i <- matrix(Dvec, nrow = num_y, ncol = nreg)

      if (njointfac != 0) {
        idioFac <- Bvec[-(1:(num_y * njointfac))]
        if (type == "countryidio_nomu") {
          idioFac <- c(idioFac, 0)
        }
        # idioPos <- TRUE
        jointFac <- matrix(Bvec[1:(num_y * njointfac)], ncol = njointfac)
        # neu
        # Unterscheidung i =1 und sonst
        # jointFac <- matrix(0, nrow = num_y, ncol = njointfac)
        # jointFac_select <- lower.tri(jointFac, diag = T)
        # jointFac[jointFac_select] <- Bvec[1:sum(jointFac_select)]
        # idioFac <- Bvec[-(1:sum(jointFac_select))]
      } else {
        idioFac <- Bvec # muss fuer regs geaendert werden bzw. Bvec oben neu definieren (D abziehen)
      }
      # idioPos <- all(idioFac > 0)
      # gewekePos <- TRUE

      # if(i == 1 & njointfac != 0){
      #   gewekeMatrix <- jointFac[1:njointfac, 1:njointfac] # muss verallgemeinert werden fuer num_y < njointfac
      #   gewekePos <- all(diag(gewekeMatrix, njointfac) > 0) # muss verallgemeinert werden fuer num_y < njointfac
      #   #gewekePos <- TRUE
      # }
      # if(gewekePos & idioPos){valid <- TRUE}
      # if(!identification){valid <- TRUE}
      # valid <- TRUE
      # ident_control <- ident_control + 1
      # }

      # if(ident_block){break}
      if (njointfac != 0) {
        Bjoint[((i - 1) * num_y + 1):(i * num_y), ] <- jointFac
      }
      Bidio[, i] <- idioFac


      if (nreg != 0) {
        DSTORE[(1 + num_y * (i - 1)):(i * num_y), (1 + nreg * (i - 1)):(i * nreg), iter] <- D_i
      }
    }
    # browser()


    # if (!ident_block) {
    if (type == "allidio") {
      if (njointfac != 0) {
        B <- cbind(Bjoint, diag(c(Bidio)))
      } else {
        B <- diag(c(Bidio))
      }
    } else if (type == "countryidio" | type == "countryidio_nomu") {
      BidioMat <- as.matrix(Matrix::bdiag(plyr::alply(array(Bidio, c(num_y, 1, NN)), 3)))
      if (njointfac != 0) {
        B <- cbind(Bjoint, BidioMat)
      } else {
        B <- BidioMat
      }
    }


    BSTORE[, , iter] <- B

    if (nreg != 0) {
      D <- DSTORE[, , iter]
    } else {
      D <- matrix(0, nrow = N_num_y)
      wReg <- matrix(0, ncol = TT)
    }


    #### GIBBS PART: Sampling VCOV matrix or Adjustment-Matrix A
    u <- compute_residuals(yObs, fPost, B, D, wReg)
    # uSTORE[,, iter] <- u
    if (VdiagEst) {
      tmp_V <- sample_V(VhatDiagScale, VhatArrayBdiagByTimeFix, VhatFix, u,
                        TT, NN_TT, num_y, availableObs_crossSection,
                        alpha0, beta0)
      VhatSqrt       <- tmp_V$V_hat_sqrt
      VSTORE[, iter] <- tmp_V$V_diag
    }

    if (sampleA) {
      uSplit <- lapply(split(u, matrix(rep(1:NN, each = TT_num_y), ncol = TT, byrow = T)), matrix, ncol = TT)

      if (countryA) {
        utu_country <- utuSum(uSplit = uSplit, VhatSqrt = VhatSqrt, N = NN, TT = TT, npara = num_y)$sumUtu_individ
        if (diagA) {
          uSum_country <- lapply(utu_country, diag)
          Adiag_country <- lapply(1:NN, \(xx)  sapply(uSum_country[[xx]], \(x) LaplacesDemon::rinvgamma(n = 1, shape = shape0 + 0.5 * sum(availableObs_crossSection[xx, ]), scale = rate0 + 0.5 * x)))
          A_countryArray <- array(unlist(lapply(Adiag_country, diag)), c(num_y, num_y, NN))
        } else {
          A_country <- lapply(1:NN, \(x) LaplacesDemon::rinvwishart(nu = sum(availableObs_crossSection[x, ]) + nu0, S = Psi0 + utu_country[[x]]))
          A_countryArray <- array(unlist(A_country), c(num_y, num_y, NN))
        }
        ASTORE[, , iter, ] <- A_countryArray
      } else {
        utu <- utuSum(uSplit = uSplit, VhatSqrt = VhatSqrt, N = NN, TT = TT, npara = num_y)$sumUtu_total
        utu <- (utu + t(utu)) / 2
        if (diagA) {
          # uSplitSqrt <- lapply(uSplit,\(x) x^2)
          # uSum <- apply(sapply(uSplitSqrt, \(x) apply(x,1,sum)),1,sum)
          uSum <- diag(utu)
          A_diag <- sapply(uSum, \(x) invgamma::rinvgamma(n = 1, shape = shape0 + 0.5 * (NN_TT_avail), rate = rate0 + 0.5 * x))
          A <- diag(A_diag)
        } else if (scaleA) {
          uSum <- sum(diag(utu))
          A_scale <- invgamma::rinvgamma(n = 1, shape = shape0 + 0.5 * (num_y * NN_TT_avail), rate = rate0 + 0.5 * uSum)
          A <- diag(rep(A_scale, num_y))
        } else {
          A <- LaplacesDemon::rinvwishart(nu = NN_TT_avail + nu0, S = Psi0 + utu)
        }
        ASTORE[, , iter] <- A
      }
    }
    # }
    ############################################################################
    ## GIBBS sampler Iteration ENDE
    ############################################################################
    # info <- sprintf("%d%% done", round((iter/itermax)*100))
    # setWinProgressBar(pb, iter/itermax*100, label=info)
    ############################################################################
    ## GIBBS sampler SPEICHERN
    ############################################################################
    # if(iter %% 10000 == 0){print(iter)}
    if ((iter == 100 & storePath != "none") | (iter == itermax & storePath != "none") | (storePath != "none" & iter %% storeUnit == 0)) {
      store_count <- store_count + 1
      if (store_count == 1) {
        dir.create(storePath_adj, recursive = T)
      }
      if (VdiagEst) {
        saveRDS(list(f = fSTORE, B = BSTORE, D = DSTORE, V = VSTORE, blockCount = block_count, errorMsg = msg_error_kf, initials = initials),
          file = paste0(storePath_adj, "/", "pj", njointfac, "_B", round(initials$B0[1, 1, 1], 2), "_Om", initials$Omega0[1, 1], "_D", initials$D[1, 1, 1], "_OmD", Omega0[dim(Omega0)[1], dim(Omega0)[1]], "_V", Vstart, "_alpha", initials$alpha0, "_beta", initials$beta0, "_IO", incObsNew, ".rds")
          # STORE WITH uSTORE:
          # saveRDS(list(f = fSTORE, B = BSTORE, D = DSTORE, V = VSTORE, u = uSTORE, blockCount = block_count,  errorMsg = errorMsg, initials = initials)
          #         , file = paste0(storePath_adj,"/", "pj",njointfac,"_B",round(initials$B0[1,1,1],2),"_Om",initials$Omega0[1,1],"_D",initials$D[1,1,1],"_OmD",Omega0[dim(Omega0)[1], dim(Omega0)[1]],"_V",Vstart, "_alpha",initials$alpha0,"_beta",initials$beta0,"_IO",incObsNew,".rds") )
        )
      } else if (sampleA) {
        saveRDS(list(f = fSTORE, B = BSTORE, D = DSTORE, A = ASTORE, blockCount = block_count, errorMsg = msg_error_kf, initials = initials),
          file = paste0(storePath_adj, "/", "cs", covScale, "_pj", njointfac, "_B", round(initials$B0[1, 1, 1], 2), "_Om", initials$Omega0[1, 1], "_D", initials$D0[1, 1, 1], "_OmD", Omega0[dim(Omega0)[1], dim(Omega0)[1]], "_A", initials$A[1, 1], "_Psi", initials$Psi0[1, 1], "_nu", initials$nu0, "_IO", incObsNew, ".rds")
        )
      } else {
        saveRDS(list(f = fSTORE, B = BSTORE, D = DSTORE, A = ASTORE, blockCount = block_count, errorMsg = msg_error_kf, initials = initials),
          file = paste0(storePath_adj, "/", "cs", covScale, "_pj", njointfac, "_B", round(initials$B0[1, 1, 1], 2), "_Om", initials$Omega0[1, 1], "_D", initials$D0[1, 1, 1], "_OmD", Omega0[dim(Omega0)[1], dim(Omega0)[1]], "_IO", incObsNew, ".rds")
        )
      }
    }

    print(iter)
  }

  # close(pb)
  return(list(f = fSTORE, B = BSTORE, D = DSTORE, A = ASTORE, V = VSTORE, blockCount = block_count, errorMsg = msg_error_kf, initials = initials))
  # RETURN WITH uSTOREÖ
  # return(list(f = fSTORE, B = BSTORE, D = DSTORE, A = ASTORE, V = VSTORE, u = uSTORE, blockCount = blockCount, errorMsg = errorMsg, initials = initials))
}
