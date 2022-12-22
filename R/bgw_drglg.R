#' @title bgw_drglg
#'
#' @description An r wrapper to call drglg_c, which in turn is a wrapper
#'   to call the Fortran subroutine drglg.
#'   drglg is the BGW iteration driver for performing statistical parameter estimation
#'
# Documentation directly from drglg:
# ! D....... SCALE VECTOR.
# ! DR...... DERIVATIVES OF R AT X.
# ! IV...... INTEGER VALUES ARRAY.
# ! LIV..... LENGTH OF IV... LIV MUST BE AT LEAST P + 90.
# ! LV...... LENGTH OF V...  LV  MUST BE AT LEAST
# !              105 + P*(2*P+16) + 2*N + 4*PS.
# ! N....... TOTAL NUMBER OF RESIDUALS.
# ! ND...... LEADING DIMENSION OF DR -- MUST BE AT LEAST PS.
# ! NN...... LEAD DIMENSION OF R, RD.
# ! P....... NUMBER OF PARAMETERS (COMPONENTS OF X) BEING ESTIMATED.
# ! PS...... NUMBER OF NON-NUISANCE PARAMETERS.
# ! R....... RESIDUALS (OR MEANS -- FUNCTIONS OF X) WHEN  DRGLG IS CALLED
# !          WITH IV(1) = 1.
# ! RD...... RD(I) = HALF * (G(I)**T * H(I)**-1 * G(I)) ON OUTPUT WHEN
# !          IV(RDREQ) IS 2, 3, 5, OR 6.   DRGLG SETS IV(REGD) = 1 IF RD
# !          IS SUCCESSFULLY COMPUTED, TO 0 IF NO ATTEMPT WAS MADE
# !          TO COMPUTE IT, AND TO -1 IF H (THE FINITE-DIFFERENCE HESSIAN)
# !          WAS INDEFINITE.  BEFORE CONVERGENCE, RD IS ALSO USED AS
# !          TEMPORARY STORAGE.
# ! RHO..... COMPUTES INFO ABOUT OBJECTIVE FUNCTION.
# ! RHOI.... PASSED WITHOUT CHANGE TO RHO.
# ! RHOR.... PASSED WITHOUT CHANGE TO RHO.
# ! V....... FLOATING-POINT VALUES ARRAY.
# ! X....... PARAMETER VECTOR BEING ESTIMATED (INPUT = INITIAL GUESS, OUTPUT = BEST VALUE FOUND).
#' @param d     Scaling vector
#' @param dr    Derivative of the choice probability model wrt x
#' @param iv    BGW internal vector of integer values
#' @param liv   Length of iv.
#' @param lv    Length of v.
#' @param n     Dimension of vector (r) of generalized residuals for the model
#' @param nd    Leading dimension dr. Must be at least ps.
#' @param nn    Leading dimension of r, rd
#' @param p     Dimension of x (as well as d, g) = number of parameters being estimated
#' @param ps    Number of non-nuisance parameters (= p in this implementation)
#' @param r     Vector of generalized residuals for the model
#' @param rd    Vector of storage space for regression diagonostics (not currently used)
#' @param v     BGW internal vector of numeric values
#' @param rhoi  Vector of integers for use by user (not currently used)
#' @param rhor  Vector of numeric values for use by user (not currently used)
#' @param x     Parameter vector for which the objective function is being minimized
#' @param i_itsum     Variable for passing itsum instruction back to bgw_mle
#'
#' @return out List of return values.
# d_vec    <- out[[1]]
# dr_vec   <- out[[2]]
# iv       <- out[[3]]
# r_vec    <- out[[4]]
# rd_vec   <- out[[5]]
# v        <- out[[6]]
# x        <- out[[7]]
# rhoi     <- out[[8]]
# rhor     <- out[[9]]
bgw_drglg <- function(d, dr, iv, liv, lv, n, nd, nn, p, ps, r, rd,
                         v, x, rhoi, rhor, i_itsum) {
    .Call(drglg_c, as.double(d), as.double(dr), as.integer(iv), as.integer(liv),
          as.integer(lv), as.integer(n),
          as.integer(nd), as.integer(nn), as.integer(p), as.integer(ps),
          as.double(r), as.double(rd), as.double(v), as.double(x),
          as.integer(rhoi), as.double(rhor), as.integer(i_itsum))
}
