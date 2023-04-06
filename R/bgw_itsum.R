#' @title bgw_itsum
#'
#' @description Prints iteration summary, info on initial and final x.
#'
#' @param d  Scaling vector
#' @param g  Gradient of objective function (negative-log-likelihood) wrt x
#' @param iv BGW internal vector of integer values
#' @param v  BGW internal vector of numeric values
#' @param x  Parameter vector for which the objective function is being minimized
#' @param p  Dimension of x
#' @param betaIsNamed Logical.
#' @param betaNames Character vector. If available, has beta parameter names.
#' @param i_itsum Code from caller to select specific options
#'
#' @return iv There are iv values changed inside bgw_itsum. The iv vector is the integer workspace for Fortran BGW.
##------------------------------------------------------------------------
# November 1, 2022
# The initial version(s) of bgw_itsum were focused on replicating (to the
# degree possible) the same functionality as ditsum in Fortran BGW.
##------------------------------------------------------------------------
bgw_itsum <- function(d, g, iv, v, x, p, betaIsNamed, betaNames=NULL, i_itsum) {
  covmat = 26
  covreq = 15
  outlev = 19
  niter = 31
  prntit = 39
  nfcall = 6
  ngcall = 30
  ngcov = 53
  nfcov = 52
  needhd = 36
  f0 = 13
  f = 10
  fdif = 11
  fdh = 74
  h = 56
  hc = 71
  preduc = 7
  prunit = 21
  nreduc = 6
  reldx = 17
  dstnrm = 2
  sused = 64
  rcond = 53
  x0prt = 24
  stppar = 5
  statpr = 23
  solprt = 22
  ##------------------------------------------------------------------------
  modelS = list()
  modelS[[1]] = '  G  '
  modelS[[2]] = '  S  '
  modelS[[3]] = ' G-S '
  modelS[[4]] = ' S-G '
  modelS[[5]] = 'G-S-G'
  modelS[[6]] = 'S-G-S'
  ##------------------------------------------------------------------------
  ## This perhaps should be moved to a function somewhere
  ## Convergence codes and messages
  bgw_message_table = list()
  bgw_message_table[[1]]  = "       ***** Evaluate function ***** \n"
  bgw_message_table[[2]]  = "       ***** Evaluate gradient ***** \n"
  # iv[3] to iv[6] are "good convergence"
  bgw_message_table[[3]]  = "       ***** X-convergence ***** \n"
  bgw_message_table[[4]]  = "       ***** Relative function convergence ***** \n"
  bgw_message_table[[5]]  = "       ***** X- and relative function convergence ***** \n"
  bgw_message_table[[6]]  = "       ***** Absolute function convergence ***** \n"
  # iv[7] and iv[8] are "bad convergence"
  bgw_message_table[[7]]  = "       ***** Singular convergence ***** \n"
  bgw_message_table[[8]]  = "       ***** False convergence ***** \n"
  # iv[9] to iv[11] are "premature stopping of search"
  bgw_message_table[[9]]  = "       ***** Function evaluation limit ***** \n"
  bgw_message_table[[10]] = "       ***** Iteration limit ***** \n"
  bgw_message_table[[11]] = "       ***** stopx ***** \n"
  # Note:  The following occurs on entry to ditsum if iv[1] > 62
  # This manifests as iv1 >= 12, means that a problem occurred at the START.
  #
  # If iv[1] > 61, iv1 = iv[1] - 51
  # iv[1] = 63 => iv1 = 12
  bgw_message_table[[12]] = "       ***** Initial f(x) cannot be computed ***** \n"
  # iv[1] = 64 => iv1 = 13
  bgw_message_table[[13]] = "       ***** Bad parameters to assess ***** \n"
  # iv[1] = 65 => iv1 = 14
  bgw_message_table[[14]] = "       ***** Gradient could not be computed ***** \n"
  # iv[1] = 66 => iv1 = 15
  bgw_message_table[[15]] = "       ***** Inconsistent dimensions (Could be N, ND, P, or LIV. )***** \n"
  # iv1 >= 16 (i.e., iv[1] >= 67) simply writes out iv[1]
  # These values would typically occur when the user overwrites
  # a default value with an illegal value. We will fill in the details.
  # iv[1] = 67 => iv1 = 16
  # This can be set two different ways in dparck. The following should cover it.
  bgw_message_table[[16]] = " ***** Internal error specifying alg in divset and/or dparck. It should be 1. ***** \n"

  bgw_success = list()
  bgw_success[[1]] = 0
  bgw_success[[2]] = 0
  bgw_success[[3]] = 0
  bgw_success[[4]] = 0
  bgw_success[[5]] = 0
  bgw_success[[6]] = 0
  bgw_success[[7]] = 1
  bgw_success[[8]] = 1
  bgw_success[[9]] = 1
  bgw_success[[10]] = 1
  bgw_success[[11]] = 1
  bgw_success[[12]] = 1
  bgw_success[[13]] = 1
  bgw_success[[14]] = 1
  bgw_success[[15]] = 1
  # In DITSUM IV(1) < 2 or > 15	Write out value of IV(1)
  ##------------------------------------------------------------------------
  # if (!betaIsNamed) betaNames <- paste0("c", 1:p)

  goto120 = FALSE
  goto390 = FALSE
  goto430 = FALSE
  goto460 = FALSE
  #
  # The following is a placeholder to indicate what
  # happens in Fortran DITSUM
  # if (iv[prunit] == 0) return(iv)
  #
  iv1 <- iv[1]
  ol <- iv[outlev] # outlev = 19
  if (iv1 > 62) iv1 <- iv1 - 51
  if ((iv1 < 2)|(iv1 > 15)) {
    # We need to document and characterize this in greater detail.
    print(" Possible internal error...")
    print(paste(" ***** iv[1] = ", iv[1]))
    return(iv)
  }

  if (iv1 < 12) {
    # Occurs: During cases after a successful start.
    # 1. During standard iterations (iv[1]-iv[2]).
    # 2. If some type of convergence has occurred (iv[3]-iv[8])).
    # 3. If the search is being stopped (iv[9]-iv[11]).
    #
    if ( (iv1 == 2)&&(iv[niter] == 0) ) goto390 <- TRUE
    if ((ol == 0) && !goto390 ) goto120 = TRUE
    if ((iv1 >=10)&&(iv[prntit]==0) && !goto390 ) goto120 = TRUE
    if ( (!goto390) && (!goto120) ) {
      if (iv1 <= 2) {
        iv[prntit] <- iv[prntit] + 1
        if (iv[prntit] < abs(ol)) return(iv)
      }
      # Line 10
      nf <- iv[nfcall] - abs(iv[nfcov])
      iv[prntit] <- 0
      reldf <- 0.
      preldf <- 0.
      oldf <- max( abs(v[f0]), abs(v[f]) )
      if (oldf > 0.) {
        reldf <- v[fdif]/oldf
        preldf <- v[preduc]/oldf
      }
      # if ( iv[needhd] == 1) {
      #   if (ol > 0) print("    it    nf      F          RELDF    PRELDF    RELDX   MODEL  STPPAR   D*step   NPRELDF")
      #   if (ol < 0) print("    it    nf      F          RELDF    PRELDF    RELDX   MODEL  STPPAR")
      # }
      # iv[needhd] <- 0
      nreldf <- 0.
      if (oldf > 0.) nreldf <- v[nreduc]/oldf
      m <- iv[sused]
      # WRITE(pu,100) stuff
      # DSB NOTE: We have not figured out why the last iteration information is printed TWICE!
      # print(iv1)
      if (iv[1] == 2) {
        if (ol > 0) cat(format(iv[niter],nsmall=0,trim=FALSE,width=6), format(nf,nsmall=0,trim=FALSE,width=5),
          formatC(v[f],format="e",digits=9),
          formatC(reldf,format="e",digits=3),
          formatC(preldf,format="e",digits=3),
          formatC(v[reldx],format="e",digits=2),
          modelS[[m]],
          formatC(v[stppar],format="e",digits=2),
          formatC(v[dstnrm],format="e",digits=2),
          formatC(nreldf,format="e",digits=3))
        if (ol < 0) cat(format(iv[niter],nsmall=0,trim=FALSE,width=6), format(nf,nsmall=0,trim=FALSE,width=5),
                        formatC(v[f],format="e",digits=9),
                        formatC(reldf,format="e",digits=3),
                        formatC(preldf,format="e",digits=3),
                        formatC(v[reldx],format="e",digits=2),
                        modelS[[m]],
                        formatC(v[stppar],format="e",digits=2))
      # Original FORMAT(i6,i5,e10.3,2E9.2,e8.1,a3,a4,2E8.1,e9.2)
      # FORMAT(i6,i5,e15.8,2E9.2,e8.1,a3,a4,2E8.1,e9.2)
      cat("\n")
      }
      goto120 <- TRUE
    }
  } else {
    goto120 <- TRUE
  }
  #
  if (goto120) {
    if (iv1 <= 2) return(iv)
    i <- iv[statpr]
    if ( (i == -1)|(i + iv1 < 0) ) {
      goto460 =TRUE
    } else {
      cat("\n")
      if (iv1 == 1) return(iv)
      if (iv1 == 2) return(iv)
      if ( (iv1 >=3)&&(iv1 <=15) ) {
        cat(bgw_message_table[[iv1]])
      }
      if ( (iv1 >=3)&&(iv1 <=11) ) goto430 = TRUE
      if (iv1 == 12) goto390 = TRUE
      if ( (iv1 == 13)|(iv1 == 15) ) return(iv)
      if (iv1 ==14) {
        if (iv(niter) > 0) { goto460 = TRUE
        } else { goto390 = TRUE
      }
      }
    }
  }
  #
  if (goto390) {
    if ( iv[x0prt] != 0 ) {
      cat("\n")
      # df <- data.frame('XName' = betaNames,
      #                  'Initial_X_i' = x, 'D_i' = d)
      # colnames(df) <- c('XName','Initial X(i)','D(i)')
      # writeLines(paste0("    ",capture.output(print(df))))
      # if (betaIsNamed) {
      #     df <- data.frame('XName' = betaNames,
      #                  'Initial_X_i' = x)
      #     colnames(df) <- c('XName','Initial X(i)')
      # } else {
      #   df <- data.frame('Initial_X_i' = x)
      #   colnames(df) <- c('Initial X(i)')
      # }
      if (betaIsNamed) {
        df <- data.frame('BetaName' = betaNames,'Initial_X_i' = x, 'D_i' = d)
        colnames(df) <- c('BetaName','InitialBeta(i)','D(i)')
      } else {
        df <- data.frame('Initial_X_i' = x,'D_i' = d)
        colnames(df) <- c('InitialBeta(i)','D(i)')
      }
      writeLines(paste0("    ",capture.output(print(df))))
      rm(df)
      cat("\n")
    }
    v[dstnrm] <- 0.
    v[fdif] <- 0.
    v[nreduc] <- 0.
    v[preduc] <- 0.
    v[reldx] <- 0.
    if (iv1 >= 12) return(iv)
    iv[needhd] <- 0
    iv[prntit] <- 0
    if (ol == 0) return(iv)
    if (ol > 0) cat("    it    nf     F            RELDF    PRELDF    RELDX    MODEL stppar   D*step   NPRELDF")
    if (ol < 0) cat("    it    nf     F            RELDF    PRELDF    RELDX    MODEL stppar")
    # cat("    it    nf     F            RELDF    PRELDF    RELDX    MODEL stppar   D*step   NPRELDF")
    cat("\n")
    cat(format(iv[niter],nsmall=0,trim=FALSE,width=6), format(iv[nfcall],nsmall=0,trim=FALSE,width=5),
        formatC(v[f],format="e",digits=9))
    cat("\n")
  }

  if (goto430) {
    goto460 <- TRUE
    # goto430
    if (iv[statpr] > 0) {
      oldf <- max( abs(v[f0]), abs(v[f]) )
      preldf = 0.
      nreldf = 0.
      if (oldf > 0.) {
        preldf <- v[preduc]/oldf
        npreldf <- v[nreduc]/oldf
      }
      nf <- iv[nfcall] - iv[nfcov]
      ng <- iv[ngcall] - iv[ngcov]
      # if (i_itsum == 1) {
      cat("\n")
      cat("       FUNCTION    ",formatC(v[f],format="e",digits=9),
          " RELDX       ", formatC(v[reldx],format="e",digits=3),
         "\n       PRELDF      ",formatC(preldf,format="e",digits=3),
            "       NPRELDF     ",formatC(npreldf,format="e",digits=3),
         "\n       func. evals ",format(nf,digits=8),
         "              grad. evals ",format(ng,digits=8),"\n")
      #}
    }
  }

  if (goto460) {
    # goto460
    if (iv[solprt] == 0) {
      cat("\n")
      return(iv)
    }
    # iv[covmat] = iv[26] = location of variance-covariance matrix
    if (iv[covmat] > 0) {
      cat("\n")
      if (iv[covreq] == 3)  cat('       vcHessianMethod = Gauss-Newton/BHHH \n')
      if (iv[covreq] == 2)  cat('       vcHessianMethod = finite differences (gradient evals) \n')
      if (iv[covreq] == -2) cat('       vcHessianMethod = finite differences (function evals) \n')
      cat('       Estimated upper bound on reciprocal of Euclidean condition number of vcHessian: ',v[rcond]^2,' \n')
      cat('       (Condition number = ratio of smallest singular value to largest singular value.)\n')
      cat('       Value of unit roundoff = machine epsilon for comparison purposes:               ',.Machine$double.eps,' \n')
      # cat('       Value of iv[fdh] = ',iv[fdh],'\n')
      # cat('       Value of iv[hc]  = ',iv[hc],'\n')
      # cat('       Value of iv[h]   = ',iv[h],'\n')


      vc_vec <- v[iv[26]:(iv[26] + p*(p+1)/2 - 1)]
      vc = matrix(0,p,p);
      iii <- 0 ;
      for (ii in 1:p) {
        for (jj in 1:ii) {
          iii <- iii + 1
          vc[ii,jj] = vc_vec[iii];
          vc[jj,ii] = vc_vec[iii];
        }
      }
      iv[needhd] <- 1
      cat("\n")
      se <- rep(0,p)
      tstat <- rep(0,p)
      for (i in 1:p) {
        se[i] <- sqrt(vc[i,i])
        tstat[i] <- x[i]/se[i]
      }
      if (betaIsNamed) {
         df <- data.frame('XName' = betaNames,
                       'FinalX(i)' = x, 's.e.(i)' = se, 't.rat.(0)' = tstat, 'G(i)' = g,'D(i)'=d)
         colnames(df) <- c('BetaName','FinalBeta(i)','s.e.(i)','t.rat.(0)','G(i)','D(i)')

      } else {
        df <- data.frame('FinalX(i)' = x, 's.e.(i)' = se, 't.rat.(0)' = tstat, 'G(i)' = g,'D(i)'=d)
        colnames(df) <- c('FinalBeta(i)','s.e.(i)','t.rat.(0)','G(i)','D(i)')
      }
      writeLines(paste0("    ",capture.output(print(df))))
      rm(se,tstat,df)
    } else {
      iv[needhd] <- 1
      # df <- data.frame('XName' = betaNames,
      #                  'FinalX(i)' = x, 'D(i)' = d, 'G(i)' = g)
      # colnames(df) <- c('XName','FinalX(i)','D(i)','G(i)')
      # writeLines(paste0("    ",capture.output(print(df))))
      # rm(df)
      rdreq <- 57
      if (iv[rdreq] == 0) {
        cat('       No variance-covariance matrix computation was requested. \n')
        cat("\n")
      } else {
        cat("\n")
        if (iv[covreq] == 3) cat('       Requested vcHessianMethod = Gauss-Newton/BHHH \n')
        if (iv[covreq] == 2) cat('       Requested vcHessianMethod = finite differences (gradient evals) \n')
        if (iv[covreq] == -2) cat('      Requested vcHessianMethod = finite differences (function evals) \n')
        if ( iv[covmat] == -1 ) {
          cat('       ***** Unable to compute a variance-covariance matrix... ***** \n')
          cat('       ***** Requested vcHessian matrix was indefinite         *****  \n')
        }
        if ( iv[covmat] == -2 ) {
          cat('      ***** Unable to compute a variance-covariance matrix....                                          ***** \n')
          cat('      ***** Too many beta values were rejected during finite-difference computation of vcHessian matrix *****  \n')
        }
        cat('       Estimated upper bound on reciprocal of Euclidean condition number of vcHessian: ',v[rcond]^2,' \n')
        cat('       (Condition number = ratio of smallest singular value to largest singular value.)\n')
        cat('       Value of unit roundoff = machine epsilon for comparison purposes:               ',.Machine$double.eps,' \n')
        # cat('       Value of iv[fdh] = ',iv[fdh],'\n')
        # cat('       Value of iv[hc]  = ',iv[hc],'\n')
        # cat('       Value of iv[h]   = ',iv[h],'\n')
        cat("\n")
      }

      if (betaIsNamed) {
         df <- data.frame('XName' = betaNames,
                       'FinalX(i)' = x, 'G(i)' = g,'D(i)'=d)
         colnames(df) <- c('BetaName','FinalBeta(i)','G(i)','D(i)')
      } else {
        df <- data.frame('FinalX(i)' = x, 'G(i)' = g,'D(i)'=d)
        colnames(df) <- c('FinalBeta(i)','G(i)','D(i)')
      }
      writeLines(paste0("    ",capture.output(print(df))))
      rm(df)
      # tmp <- c(x,d,g)
      # # tmp <- matrix(tmp, nrow = length(tmp), ncol = 3, dimnames=list(betaNames,c('Final X(i)','D(i)','G(i)')))
      # tmp <- matrix(tmp,nrow = length(tmp), ncol = 3)
      # print(tmp)
      # rm(tmp)
      # output = cbind(c(1:p),x,d,g)
      # colnames(output) <- c('i','Final X(i)','D(i)','G(i)')
      # rownames(output) <- betaNames
      # apollo_print(output) #print(output, digits = 6)
    }
    cat("\n")
  }
  return(iv)
}
