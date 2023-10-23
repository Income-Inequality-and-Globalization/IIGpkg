#' Title
#'
#' @param itermax Number of Gibbs iterations
#' @param identmax Nicht mehr noetig (war fuer altes Sampling proceduere - so lange samplen bis Parameterrestriktionen eingehalten)
#' @param npara Anzahl der Parameter (bei uns immer 3)
#' @param nreg Anzahl der  Regressoren (1 bis 3)
#' @param njointfac Anzahl der gemeinsamen Faktoren (zuletzt immer 0)
#' @param yObs Daten (a, q, mu)
#' @param c_ Konstante der Transition-equation (keine Konstante im Modell - daher nicht relevant)
#' @param D Startwert fuer die Koeffizienten der Regressorn (stacked: npara*N x nreg*N), N = Anzahl Laender
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
#' @param VhatDiagScale (wahrscheinlich irrelevant) = TRUE, wenn fuer jedes Diagonalelement der GMM- varianz ein Skalierungsfaktor geschaetzt werden soll. In Gibb2_SM_SA_cluster muss dann Diag_VCOV = TRUE gesetzt werden, damit die Nebendiagonalelemente der GMM-Varianz = 0 gesetzt werden. Wenn VDiagScale = TRUE sein soll, muss auch VDiagEst = TRUE sein.
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
#' @return
#' @export

GibbsSSM_2 <- function(itermax = 15000, identmax = 5000, npara, nreg, njointfac, yObs, c_, D, D0, B, Phi, Q, type = "allidio", initX, initP, initU, wRegSpec, wReg, uReg = NULL,
                       B0, Omega0, selectR, selectC, Vhat, incObsOld = 100000, incObsNew = 100000, covScale, VhatDiagScale, VhatDiagScale_start = NULL, VdiagEst, alpha0, beta0, countryA, A, scaleA, diagA, Psi0, nu0, shape0, rate0, storePath = "none", fPost, sampleA, identification = T, storeUnit = 10000) {
  # Vhat muss (npara, npara, N * T) Array sein. Start: Alle Zeitpunkte fuer das erste Land.
  # Funktioniert noch nicht fuer npara < njointfac (Geweke)
  # Modell ohne Konstante
  initials <- c(as.list(environment()), list()) # speichert die Initialisierung

  storeCount <- 0 # zaehlt wie oft bereits zwischengespeichert wurde
  Vhat <- incObsOld / incObsNew * Vhat

  if (missing(fPost)) {
    Simf <- TRUE # Simf = TRUE, wenn FFBS durchgefuerht werden soll (immer der Fall, wenn fPost bei der Initialisierung fehlt)
  } else {
    Simf <- FALSE
  }

  if (missing(wRegSpec)) {
    wRegSpec <- 0
  }

  errorMsg <- NULL # Zum Speichern eines Fehlers im KF

  nfac <- dim(B)[2] # number of factors
  Nnpara <- dim(B)[1] # N * npara

  N <- Nnpara / npara # Anzahl der Laender
  TT <- dim(yObs)[2] # Anzahl der Zeitpunkte
  fSTORE <- array(0, dim = c(nfac, TT, itermax)) # speichert den Output von FFBS
  # cSTORE <- matrix(rep(0,Nnpara*itermax), ncol = itermax)
  BSTORE <- array(0, dim = c(Nnpara, nfac, itermax)) # speichert die Gibbs-Zuege der Ladungen auf die latenten Faktoren
  # uSTORE <- array(0, dim = c(Nnpara, TT, itermax) )
  if (nreg != 0) {
    DSTORE <- array(0, dim = c(Nnpara, N * nreg, itermax)) # speichert die Koeffizienten der Makroregressor-Koeffizienten
    dimOmega0 <- dim(Omega0)[1]
    OmegaD0 <- Omega0[dimOmega0, dimOmega0] # speichert Startwert der Priorvarianz fuer die Makroregressor-Koeffizienten. Da ich Omega0 immer diagonal waehle und fuer jeden Makroregressor-Koeffizienten die gleiche Varianz waehle, reicht es hier einen Parameter zu speichern.
  } else {
    DSTORE <- NULL
    D0 <- NULL
    OmegaD0 <- NULL
  }

  if (countryA) {
    ASTORE <- array(0, dim = c(npara, npara, itermax, N)) # speichert die Gibbs-Zuege der Adjustmentmatrix
  } else {
    ASTORE <- array(0, dim = c(npara, npara, itermax))
  }
  VSTORE <- matrix(0, nrow = npara * N, ncol = itermax) # speichert Gibbs-Zuege für die Messfehlervarianzen, wenn diese geschaetzt werden (VDiagEst = TRUE) (wahrscheinlich fuer uns irrelevant)

  if (VdiagEst) {
    Vstart <- initials$Vhat[1, 1, 1]
  }

  if (VhatDiagScale) {
    VhatFix <- Vhat # speichert Vhat Initialisierung
    VhatArrayBdiagByTimeFix <- bdiagByTime(VhatFix, npara = npara, N = N, TT = TT, Nnpara = Nnpara) # Sortiert Vhat nach der Zeit: wie Vhat eine npara x npara x N*TT Array, aber die Elemente [1:npara, 1:npara, 1:N] sind nun die Messfehlervarianzen aller Laender zum ersten Zeitpunkt, [1:npara, 1:npara, (N+1):(2N)] alle Laender zum zweiten Zeitpunkt, usw.

    if (!missing(VhatDiagScale_start) & !is.null(VhatDiagScale_start)) {
      Vhat <- array(sapply(1:(N * TT), \(x) VhatDiagScale_start[, , x] %*% Vhat[, , x]), dim = c(npara, npara, N * TT)) # setzt Vhat im Fall VhatDiagScale = TRUE
      Vstart <- VhatDiagScale_start[1, 1, 1]
    } else {
      Vstart <- 1
    }
  }
  if (countryA) {
    A_country <- plyr::alply(replicate(N, A), 3) # Liste: enthaelt fuer jedes Land eine eigene A Matrix nach dem Vorbild des uebergebenen A's.
    A_countryArray <- array(unlist(A_country), c(npara, npara, N)) # macht aus der Liste ein Array
  }



  VhatSqrt <- array(apply(Vhat, 3, matSqrt), c(npara, npara, N * TT)) # Array mit Matrixwurzeln der Vhat-Matrizen

  NTT_available <- sum(!is.na(yObs)) / npara # N *TT ohne die fehlenden Beobachtungen
  availableObs_crossSection <- t(apply(yObs[seq(1, N * npara, 3), ], 1, \(x) !is.na(x))) # Anzahl der Beobachtungen (ohne die NAs) fuer jedes Land


  # storePath Anpassung
  if (storePath != "none" & !missing(storePath)) {
    if (VdiagEst) {
      storePath_adj <- paste0(storePath, "/", "cs", covScale, "_pj", njointfac, "_Reg", wRegSpec, "_B", round(initials$B0[1, 1, 1], 2), "_Om", initials$Omega0[1, 1], "_D", initials$D0[1, 1, 1], "_OmD", OmegaD0, "_A", initials$A[1, 1], "_Psi", initials$Psi0[1, 1], "_nu", initials$nu0, "_IO", incObsNew)
    } else if (sampleA) {
      storePath_adj <- paste0(storePath, "/", "cs", covScale, "_pj", njointfac, "_Reg", wRegSpec, "_B", round(initials$B0[1, 1, 1], 2), "_Om", initials$Omega0[1, 1], "_D", initials$D0[1, 1, 1], "_OmD", OmegaD0, "_A", initials$A[1, 1], "_Psi", initials$Psi0[1, 1], "_nu", initials$nu0, "_IO", incObsNew)
    } else {
      storePath_adj <- paste0(storePath, "/", "cs", covScale, "_pj", njointfac, "_Reg", wRegSpec, "_B", round(initials$B0[1, 1, 1], 2), "_Om", initials$Omega0[1, 1], "_D", initials$D0[1, 1, 1], "_OmD", OmegaD0, "_IO", incObsNew)
    }
  }

  # pb <- winProgressBar(title="Progress", label="0% done", min=0, max=100, initial=0)
  blockCount <- 0 # wird nicht mehr benoetigt, nicht auskommentiert, da es unten auch noch drin steht. War vom alten Sampler, um die identifizierenden Restriktionen einzuhalten.
  ident_block <- FALSE # wird nicht mehr benoetigt, nicht auskommentiert, da es unten auch noch drin steht. War vom alten Sampler, um die identifizierenden Restriktionen einzuhalten.
  iter <- 1 # aktuelle Gibbs-Iteration

  
  ## Gibbs Sampler Iterations
  while (iter <= itermax) {
    if (!ident_block) {
      
      # Erweiterung das Vhat Arrays um die Adjustmentmatrix A
      if (countryA) {
        VhatArray_A <- array(0, dim = c(npara, npara, N * TT))

        for (i in 1:N) {
          VhatArray_A[, , (TT * (i - 1) + 1):(i * TT)] <- array(apply(VhatSqrt[, , (TT * (i - 1) + 1):(i * TT)], 3, \(x) x %*% A_countryArray[, , i] %*% t(x)), c(npara, npara, TT))
        }
      } else {
        VhatArray_A <- array(apply(VhatSqrt, 3, function(X) {
          X %*% A %*% t(X)
        }), c(npara, npara, N * TT))
      }
      
      # Vhat Array (mit Adjustmentmatrix A) bzgl. der Zeit sortiert
      VhatArrayBdiagByTime <- bdiagByTime(VhatArray_A = VhatArray_A, npara = npara, N = N, TT = TT, Nnpara = Nnpara)
    
    }

    
    #### GIBBS PART: Sampling of the latent factors (FFBS)
    if (Simf) {
      if (!ident_block) {
        
        # Kalman-Filer
        invisible(capture.output(KF <- tryCatch(
          {
            RcppSMCkalman::kfMFPD(
              yObs = yObs, uReg = uReg, wReg = wReg,
              dimX = nfac, dimY = Nnpara, TT = TT,
              A = Phi, B = NULL, C = B, D = D,
              Q = Q, R = VhatArrayBdiagByTime, x00 = initX,
              u00 = initU, P00 = initP, PDSTORE = F
            )
          },
          error = function(e) {
            e$message
          }
        )))
        # try({KF <-RcppSMCkalman::kfMFPD(yObs = yObs, uReg = uReg, wReg = wReg,
        #                                                      dimX = nfac, dimY = Nnpara, TT = TT,
        #                                                      A = Phi, B = NULL, C = B, D = c_,
        #                                                      Q = Q, R = VhatArrayBdiagByTime, x00 = initX,
        #                                                     u00 = initU, P00 = initP, PDSTORE = F)})
        # if(inherits(KF, "try-error")){
        #   print("Fehler")
        #   break
        # }
      }
      
      # Speichert moegliche Fehlermeldung des Kalman-Filter und die Zwischenergebnisse zum Zeitpunkt des Fehlers zurueck
      if (is.character(KF)) {
        errorMsg <- KF
        if (storeCount == 0) {
          dir.create(storePath_adj, recursive = T)
          storeCount <- storeCount + 1
        }
        if (VdiagEst) {
          saveRDS(list(f = fSTORE, B = BSTORE, V = VSTORE, blockCount = blockCount, errorMsg = errorMsg, initials = initials),
            file = paste0(storePath_adj, "/", "B", round(initials$B0[1, 1, 1], 2), "_Omega", initials$Omega0[1, 1], "_V", initials$Vhat[1, 1, 1], "_alpha", initials$alpha0, "_beta", initials$beta0, "_IO", incObsNew, "_error", iter, ".rds")
          )
        } else {
          saveRDS(list(f = fSTORE, B = BSTORE, A = ASTORE, blockCount = blockCount, errorMsg = errorMsg, initials = initials),
            file = paste0(storePath_adj, "/", "B", round(initials$B0[1, 1, 1], 2), "_Omega", initials$Omega0[1, 1], "_A", initials$A[1, 1], "_Psi", initials$Psi0[1, 1], "_nu0", initials$nu0, "_IO", incObsNew, "_error", iter, ".rds")
          )
        }
        return(list(f = fSTORE, B = BSTORE, D = DSTORE, A = ASTORE, blockCount = blockCount, errorMsg = errorMsg))
      }
      filt_f <- KF$mfdEXP
      filt_P <- KF$mfdVAR

      ident_block <- FALSE # wird nicht mehr benoetigt, nicht auskommentiert, da es unten auch noch drin steht. War vom alten Sampler, um die identifizierenden Restriktionen einzuhalten.
      
      # Backward Sampling
      fPost <- GibbsSSM_f(TT = TT, nfac = nfac, Phi = Phi, Q = Q, filt_f = filt_f, filt_P = filt_P)
    }


    fSTORE[, , iter] <- fPost
    

    ##### GIBBS PART: Sampling of loadings (on latent factors) and partial effects (on regressors)
    # Ergebnis-Matrix fuer die gemeinsamen Ladungen (wird spaeter befuellt)
    Bjoint <- matrix(rep(0, Nnpara * njointfac), ncol = njointfac)
    # Ergebnis-Matrix fuer die idiosynkratischen Ladungen (wird spaeter befuellt)
    Bidio <- matrix(rep(0, Nnpara), ncol = N)

    for (i in 1:N) {
      # Vhat Array for cross-sectional unit (country) i
      Viarray <- VhatArray_A[, , (1 + (i - 1) * TT):(i * TT)]
      # yObs for cross-sectional unit (country) i
      yiObs <- yObs[(1 + npara * (i - 1)):(npara * i), ]
      # availableObs <- which(!is.na(yObs[1 + npara * (i-1),]))
      availableObs <- which(availableObs_crossSection[i, ])

      
      ## Posterior Momente fuer die Ladungen und die partiellen Effekte
      invOmega0 <- solve(Omega0)
      invOmega1_part2 <- sumffkronV(availableObs, npara = npara, nreg = nreg, njointfac = njointfac, i = i, fPost = fPost, wReg = wReg, Viarray = Viarray, type = type)
      invOmega1 <- invOmega0 + invOmega1_part2

      beta1_mid <- sumfyV(availableObs, npara = npara, nreg = nreg, njointfac = njointfac, i = i, fPost = fPost, wReg = wReg, yiObs = yiObs, Viarray = Viarray, type = type)

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
            if (storeCount == 0) {
              dir.create(storePath_adj, recursive = T)
              storeCount <- storeCount + 1
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
          lower <- c(rep(0, npara), rep(-Inf, npara * nreg)) # muss fuer gemeinsamen faktor angepasst werden
          upper <- rep(Inf, npara + npara * nreg) # muss fuer gemeinsamen faktor angepasst werden
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
      # muss verallgemeinert werden fuer npara < njointfac (fuer uns nicht noetig)


      # Bvec <- MASS::mvrnorm(n = 1, mu = selectR %*% beta1 + selectC, Sigma = selectR %*% Omega1 %*% t(selectR))
      # Sigma <- 0.5 * (selectR %*% Omega1 %*% t(selectR)) + 0.5 * t(selectR %*% Omega1 %*% t(selectR))
      Sigma <- 0.5 * Omega1 + 0.5 * t(Omega1)
      # Bvec <- as.numeric(tmvtnorm::rtmvnorm(n = 1, mean = as.numeric(selectR %*% beta1 + selectC), sigma = Sigma,
      #                          lower = lower, upper = upper))
       
      # Bvec <- tmvnsim::tmvnsim(1,length(upper),lower = lower, upper = upper, means = as.numeric(selectR %*% beta1 + selectC), sigma = Sigma )$samp
      # Bvec <- tmvnsim::tmvnsim(1,length(upper),lower = lower, upper = upper, means = as.numeric(beta1 + selectC), sigma = Sigma )$samp
      
      # Sampling der Ladungen bzw. partiellen Effekte
      BDsamp <- tmvnsim::tmvnsim(1, length(upper), lower = lower, upper = upper, means = as.numeric(beta1 + selectC), sigma = Sigma)$samp

      Bvec <- BDsamp[1:((njointfac + 1) * npara)]
      Dvec <- BDsamp[-(1:((njointfac + 1) * npara))]
      D_i <- matrix(Dvec, nrow = npara, ncol = nreg)

      if (njointfac != 0) {
        idioFac <- Bvec[-(1:(npara * njointfac))]
        if (type == "countryidio_nomu") {
          idioFac <- c(idioFac, 0)
        }
        # idioPos <- TRUE
        jointFac <- matrix(Bvec[1:(npara * njointfac)], ncol = njointfac)
        # neu
        # Unterscheidung i =1 und sonst
        # jointFac <- matrix(0, nrow = npara, ncol = njointfac)
        # jointFac_select <- lower.tri(jointFac, diag = T)
        # jointFac[jointFac_select] <- Bvec[1:sum(jointFac_select)]
        # idioFac <- Bvec[-(1:sum(jointFac_select))]
      } else {
        idioFac <- Bvec # muss fuer regs geaendert werden bzw. Bvec oben neu definieren (D abziehen)
      }
      # idioPos <- all(idioFac > 0)
      # gewekePos <- TRUE

      # if(i == 1 & njointfac != 0){
      #   gewekeMatrix <- jointFac[1:njointfac, 1:njointfac] # muss verallgemeinert werden fuer npara < njointfac
      #   gewekePos <- all(diag(gewekeMatrix, njointfac) > 0) # muss verallgemeinert werden fuer npara < njointfac
      #   #gewekePos <- TRUE
      # }
      # if(gewekePos & idioPos){valid <- TRUE}
      # if(!identification){valid <- TRUE}
      # valid <- TRUE
      # ident_control <- ident_control + 1
      # }

      # if(ident_block){break}
      if (njointfac != 0) {
        Bjoint[((i - 1) * npara + 1):(i * npara), ] <- jointFac
      }
      Bidio[, i] <- idioFac


      if (nreg != 0) {
        DSTORE[(1 + npara * (i - 1)):(i * npara), (1 + nreg * (i - 1)):(i * nreg), iter] <- D_i
      }
    }


    if (!ident_block) {
      if (type == "allidio") {
        if (njointfac != 0) {
          B <- cbind(Bjoint, diag(c(Bidio)))
        } else {
          B <- diag(c(Bidio))
        }
      } else if (type == "countryidio" | type == "countryidio_nomu") {
        BidioMat <- as.matrix(Matrix::bdiag(plyr::alply(array(Bidio, c(npara, 1, N)), 3)))
        if (njointfac != 0) {
          B <- cbind(Bjoint, BidioMat)
        } else {
          B <- BidioMat
        }
      }


      BSTORE[, , iter] <- B

      if (nreg != 0) {
        D <- DSTORE[, , iter]
      } else{
        D <- matrix(0, nrow = Nnpara)
        wReg <- matrix(0, ncol = TT )
      }

      
      #### GIBBS PART: Sampling VCOV matrix or Adjustment-Matrix A
      u <- yObs - B %*% fPost - D %*% wReg
      # uSTORE[,, iter] <- u
      if (VdiagEst) {
        # u <- yObs - apply(fPost, 2, function(x) {
        #   B %*% x
        # })
        if (VhatDiagScale) {
          u <- sapply(1:TT, \(x) (1 / sqrt(diag(VhatArrayBdiagByTimeFix[, , x]))) * u[, x])
        }
        ssu <- apply(u, 1, \(x) sum(x^2, na.rm = TRUE))
        TT_available <- rep(apply(availableObs_crossSection, 1, sum), each = npara)
        ssu_TT <- cbind(ssu, TT_available)
        Vdiag <- apply(ssu_TT, 1, \(x) LaplacesDemon::rinvgamma(1, shape = alpha0 + x[2] / 2, scale = beta0 + x[1] / 2))
        VdiagMat <- matrix(Vdiag, ncol = npara, byrow = T)
        VdiagMat_extend <- VdiagMat[rep(1:nrow(VdiagMat), each = TT), ]
        VdiagArray <- array(apply(VdiagMat_extend, 1, diag), dim = c(npara, npara, N * TT))
        if (VhatDiagScale) {
          VdiagArray <- array(sapply(1:(N * TT), \(x) VdiagArray[, , x] %*% VhatFix[, , x]), dim = c(npara, npara, N * TT))
        }
        VhatSqrt <- array(apply(VdiagArray, 3, matSqrt), c(npara, npara, N * TT))
        VSTORE[, iter] <- Vdiag
      }

      if (sampleA) {
        # u <- yObs - apply(fPost, 2, function(x) {
        #   B %*% x
        # })
        uSplit <- lapply(split(u, matrix(rep(1:N, each = npara * TT), ncol = TT, byrow = T)), matrix, ncol = TT)

        if (countryA) {
          utu_country <- utuSum(uSplit = uSplit, VhatSqrt = VhatSqrt, N = N, TT = TT, npara = npara)$sumUtu_individ
          if (diagA) {
            uSum_country <- lapply(utu_country, diag)
            Adiag_country <- lapply(1:N, \(xx)  sapply(uSum_country[[xx]], \(x) LaplacesDemon::rinvgamma(n = 1, shape = shape0 + 0.5 * sum(availableObs_crossSection[xx, ]), scale = rate0 + 0.5 * x)))
            A_countryArray <- array(unlist(lapply(Adiag_country, diag)), c(npara, npara, N))
          } else {
            A_country <- lapply(1:N, \(x) LaplacesDemon::rinvwishart(nu = sum(availableObs_crossSection[x, ]) + nu0, S = Psi0 + utu_country[[x]]))
            A_countryArray <- array(unlist(A_country), c(npara, npara, N))
          }
          ASTORE[, , iter, ] <- A_countryArray
        } else {
          utu <- utuSum(uSplit = uSplit, VhatSqrt = VhatSqrt, N = N, TT = TT, npara = npara)$sumUtu_total
          utu <- (utu + t(utu)) / 2
          if (diagA) {
            # uSplitSqrt <- lapply(uSplit,\(x) x^2)
            # uSum <- apply(sapply(uSplitSqrt, \(x) apply(x,1,sum)),1,sum)
            uSum <- diag(utu)
            A_diag <- sapply(uSum, \(x) invgamma::rinvgamma(n = 1, shape = shape0 + 0.5 * (NTT_available), rate = rate0 + 0.5 * x))
            A <- diag(A_diag)
          } else if (scaleA) {
            uSum <- sum(diag(utu))
            A_scale <- invgamma::rinvgamma(n = 1, shape = shape0 + 0.5 * (npara * NTT_available), rate = rate0 + 0.5 * uSum)
            A <- diag(rep(A_scale, npara))
          } else {
            A <- LaplacesDemon::rinvwishart(nu = NTT_available + nu0, S = Psi0 + utu)
          }
          ASTORE[, , iter] <- A
        }
      }
    }

    # info <- sprintf("%d%% done", round((iter/itermax)*100))
    #
    # setWinProgressBar(pb, iter/itermax*100, label=info)

    # if(iter %% 10000 == 0){print(iter)}
    if ((iter == 100 & storePath != "none") | (iter == itermax & storePath != "none") | (storePath != "none" & iter %% storeUnit == 0 & !ident_block)) {
      storeCount <- storeCount + 1
      if (storeCount == 1) {
        dir.create(storePath_adj, recursive = T)
      }
      if (VdiagEst) {
        saveRDS(list(f = fSTORE, B = BSTORE, D = DSTORE, V = VSTORE, blockCount = blockCount, errorMsg = errorMsg, initials = initials),
          file = paste0(storePath_adj, "/", "pj", njointfac, "_B", round(initials$B0[1, 1, 1], 2), "_Om", initials$Omega0[1, 1], "_D", initials$D[1, 1, 1], "_OmD", OmegaD0, "_V", Vstart, "_alpha", initials$alpha0, "_beta", initials$beta0, "_IO", incObsNew, ".rds")
          # STORE WITH uSTORE:
          # saveRDS(list(f = fSTORE, B = BSTORE, D = DSTORE, V = VSTORE, u = uSTORE, blockCount = blockCount,  errorMsg = errorMsg, initials = initials)
          #         , file = paste0(storePath_adj,"/", "pj",njointfac,"_B",round(initials$B0[1,1,1],2),"_Om",initials$Omega0[1,1],"_D",initials$D[1,1,1],"_OmD",OmegaD0,"_V",Vstart, "_alpha",initials$alpha0,"_beta",initials$beta0,"_IO",incObsNew,".rds") )
        )
      } else if (sampleA) {
        saveRDS(list(f = fSTORE, B = BSTORE, D = DSTORE, A = ASTORE, blockCount = blockCount, errorMsg = errorMsg, initials = initials),
          file = paste0(storePath_adj, "/", "cs", covScale, "_pj", njointfac, "_B", round(initials$B0[1, 1, 1], 2), "_Om", initials$Omega0[1, 1], "_D", initials$D0[1, 1, 1], "_OmD", OmegaD0, "_A", initials$A[1, 1], "_Psi", initials$Psi0[1, 1], "_nu", initials$nu0, "_IO", incObsNew, ".rds")
        )
      } else {
        saveRDS(list(f = fSTORE, B = BSTORE, D = DSTORE, A = ASTORE, blockCount = blockCount, errorMsg = errorMsg, initials = initials),
          file = paste0(storePath_adj, "/", "cs", covScale, "_pj", njointfac, "_B", round(initials$B0[1, 1, 1], 2), "_Om", initials$Omega0[1, 1], "_D", initials$D0[1, 1, 1], "_OmD", OmegaD0, "_IO", incObsNew, ".rds")
        )
      }
    }

    iter <- iter + 1
  }

  # close(pb)
  return(list(f = fSTORE, B = BSTORE, D = DSTORE, A = ASTORE, V = VSTORE, blockCount = blockCount, errorMsg = errorMsg, initials = initials))
  # RETURN WITH uSTOREÖ
  # return(list(f = fSTORE, B = BSTORE, D = DSTORE, A = ASTORE, V = VSTORE, u = uSTORE, blockCount = blockCount, errorMsg = errorMsg, initials = initials))
}
