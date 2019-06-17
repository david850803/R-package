#' Trapezoidal Integration
#' @description The numerical integration of ftn from a to b,using the trapezoidal rule with n subdivisions.
#' @param ftn a function of single variable.
#' @param a lower bound of the integration. Must be finite.
#' @param b upper bound of the integration. Must be finite and assume b > a.
#' @param n number of partitions.
#' @details The integration calculates the trapezoid area of each subdivisions and sums them up.
#' @export
#' @seealso \code{integrate}
#' @references Owen Jones,Robert Maillardet,Andrew Robinson (2014).Introduction to Scientific Programming and Simulation Using R ; Second Edition
#' @examples
#' ftn1 <-  function(x) return(4 * x^3)
#' trapezoid(ftn1, 0, 1, n = 100)


trapezoid <- function(ftn, a, b, n = 100)
{
  h <- (b-a)/n
  x.vec <- seq(a, b, by = h)
  f.vec <- sapply(x.vec, ftn)
  Total <- h * (f.vec[1]/2 + sum(f.vec[2:n]) + f.vec[n+1]/2)
  return(Total)
}

#' Adapted trapezoidal integration
#' @description  The function returns a vector of length 2 whose first element is the integration and whose second element is the number of function evaluations required.
#' @param ftn a function of single variable.
#' @param a lower bound of the integration. Must be finite.
#' @param b upper bound of the integration. Must be finite and assume b > a.
#' @param tol the partition used is recursively refined until the estimate on successive partitions differs by at most tol.
#' @param trace if trace is TRUE then intermediate results are printed.
#' @details The variable tol is the tolerant error of the approximated integration and the partitions are subdivided.
#' t.recursion is the function that recursively improves the partitions such that the summation of the errors is less than the tolerant error.
#' @export
#' @seealso \code{integrate}
#' @references Owen Jones,Robert Maillardet,Andrew Robinson (2014).Introduction to Scientific Programming and Simulation Using R ; Second Edition
#' @examples
#' ftn1 <- function(x) return(4 * x^3)
#' trapezoid_adpt(ftn1, 0, 1, 0.00001, trace = FALSE)

trapezoid_adpt <- function(ftn, a, b, tol = 1e-8, trace = FALSE)
{
  c <- (a + b)/2
  fa <- ftn(a)
  fb <- ftn(b)
  fc <- ftn(c)
  width <- c(c - a,b - c)
  I.start <- 1/2 * ( width[1] * (fa + fc) + width[2] * (fc + fb)) # Trapezoidal rule
  t.out <- t.recursion(ftn, a, b, c, fa, fb, fc,I.start, tol, 1, trace)
  t.out[2] <- t.out[2] + 3
  if (trace) {
    cat("final value is", t.out[1], "in",
        t.out[2], "function evaluations\n")
  }
  return(t.out)
}

t.recursion <- function(ftn, a, b, c, fa, fb, fc, I.old, tol, level, trace)
{
  level.max <- 100
  if (level > level.max)
  {
    cat("recursion limit reached: singularity likely\n")
    return(NULL)
  }
  else
  {
    point_left <- (a + c)/2
    point_right <- (c + b)/2
    f1 <- ftn(point_left)
    f2 <- ftn(point_right)
    width_left <- c(point_left - a, c - point_left)
    width_right <- c(point_right - c, b - point_right)
    I.left <- 1/2 * ( width_left[1] * (fa + f1) + width_left[2] * (f1 + fc))  # Trapezoid's rule for left half
    I.right <- 1/2 * ( width_right[1] * (fc + f2) + width_right[2] * (f2 + fb)) # Trapezoid's rule for right half
    I.new <- I.left + I.right       # new estimate for the integral
    f.count <- 2

    if (abs(I.new - I.old) > tol) { # I.new not accurate enough
      t.left <- t.recursion(ftn, a, c, point_left, fa, fc, f1, I.left,
                            tol/2, level + 1, trace)
      t.right <- t.recursion(ftn, c, b, point_right, fc, fb, f2, I.right,
                             tol/2, level + 1, trace)
      I.new <- t.left[1] + t.right[1]
      f.count <-  f.count + t.left[2] + t.right[2];
    } else { # we have achieved the desired tolerance
      if (trace) {
        cat("integral over [", a, ", ", b, "] is ", I.new,
            " (at level ", level, ")\n", sep = "")
      }
    }

    return(c(I.new, f.count))
  }
}

#' Trapezoidal Integration with random partitions
#' @param ftn a function of single variable.
#' @param a lower bound of the integration. Must be finite.
#' @param b upper bound of the integration. Must be finite and assume b > a.
#' @param n number of partitions.
#' @details The lengths of the partitions are not equal and the partitions are randomly selected by uniform distribution[a, b].
#' The uniform distribution has density

#' f(x) = 1/(max - min)

#' for min ≤ x ≤ max.
#' @export
#' @seealso \code{integrate}
#' @references Owen Jones,Robert Maillardet,Andrew Robinson (2014).Introduction to Scientific Programming and Simulation Using R ; Second Edition
#' @examples
#' ## You can set.seed() before using this function in order
#' ## to lead to the same result.
#' set.seed(1)
#' ftn1 <- function(x) return(4 * x^3)
#' trapezoid_rp(ftn1, 0, 1, n = 100)


trapezoid_rp <- function(ftn, a, b, n)
{
  points_add <- sort(runif((n - 1), a, b))
  points <- c(a,points_add,b)
  f.vec <- sapply(points,ftn)
  width <- NULL
  area <- NULL
  for (i in 1:n)
  {
    width[i] <- points[i + 1] - points[i]
    area[i] <- 1/2 * width[i] * (sum(f.vec[i:(i + 1)]))
  }
  Total <- sum(area)
  return(Total)
}


#' Adapted trapezoidal integration with random partitions
#' @description  The function returns a vector of length 2 whose first element is the integration and whose second element is the number of function evaluations required.
#' @param ftn a function of single variable.
#' @param a lower bound of the integration. Must be finite.
#' @param b upper bound of the integration. Must be finite and assume b > a.
#' @param tol the partition used is recursively refined until the estimate on successive partitions differs by at most tol.
#' @param trace if trace is TRUE then intermediate results are printed.
#' @details The variable tol is the tolerant error of the approximated integration and the partitions are randomly selected by uniform distribution [a, b].
#' The uniform distribution has density

#' f(x) = 1/(max - min)

#' for min ≤ x ≤ max.
#' t.recursion is the function that recursively improves the partitions such that the summation of the errors is less than the tolerant error.
#' @export
#' @seealso \code{integrate}
#' @references Owen Jones,Robert Maillardet,Andrew Robinson (2014).Introduction to Scientific Programming and Simulation Using R ; Second Edition
#' @examples
#' ## You can set.seed() before using this function in order
#' ## to lead to the same result.
#' set.seed(1)
#' ftn1 <- function(x) return(4 * x^3)
#' trapezoid_rp_adpt(ftn1, 0, 1, 0.00001, trace = FALSE)

trapezoid_rp_adpt <- function(ftn, a, b, tol = 1e-8, trace = FALSE)
{
  c <- runif(1, a, b)
  fa <- ftn(a)
  fb <- ftn(b)
  fc <- ftn(c)
  width <- c(c - a,b - c)
  I.start <- 1/2 * ( width[1] * (fa + fc) + width[2] * (fc + fb)) # Trapezoidal rule
  t.out <- t.rp.recursion(ftn, a, b, c, fa, fb, fc,I.start, tol, 1, trace)
  t.out[2] <- t.out[2] + 3
  if (trace) {
    cat("final value is", t.out[1], "in",
        t.out[2], "function evaluations\n")
  }
  return(t.out)
}

t.rp.recursion <- function(ftn, a, b, c, fa, fb, fc, I.old, tol, level, trace)
{
  level.max <- 100
  if (level > level.max)
  {
    cat("recursion limit reached: singularity likely\n")
    return(NULL)
  }
  else
  {
    point_left <- runif(1,a,c)
    point_right <- runif(1,c,b)
    f1 <- ftn(point_left)
    f2 <- ftn(point_right)
    width_left <- c(point_left - a, c - point_left)
    width_right <- c(point_right - c, b - point_right)
    I.left <- 1/2 * ( width_left[1] * (fa + f1) + width_left[2] * (f1 + fc))  # Trapezoid's rule for left half
    I.right <- 1/2 * ( width_right[1] * (fc + f2) + width_right[2] * (f2 + fb)) # Trapezoid's rule for right half
    I.new <- I.left + I.right       # new estimate for the integral
    f.count <- 2

    if (abs(I.new - I.old) > tol) { # I.new not accurate enough
      t.left <- t.rp.recursion(ftn, a, c, point_left, fa, fc, f1, I.left,
                            tol/2, level + 1, trace)
      t.right <- t.rp.recursion(ftn, c, b, point_right, fc, fb, f2, I.right,
                             tol/2, level + 1, trace)
      I.new <- t.left[1] + t.right[1]
      f.count <-  f.count + t.left[2] + t.right[2];
    } else { # we have achieved the desired tolerance
      if (trace) {
        cat("integral over [", a, ", ", b, "] is ", I.new,
            " (at level ", level, ")\n", sep = "")
      }
    }

    return(c(I.new, f.count))
  }
}


#' Monte-Carlo integration using the hit and miss method
#' @param ftn a function of single variable.
#' @param a lower bound of the integration. Must be finite.
#' @param b upper bound of the integration. Must be finite and assume b > a.
#' @param f.min lower bound of the ftn over the range [a,b].
#' @param f.max upper bound of the ftn over the range [a,b].
#' @param n number of samples used in the estimation.
#' @export
#' @seealso \code{integrate}
#' @references Owen Jones,Robert Maillardet,Andrew Robinson (2014).Introduction to Scientific Programming and Simulation Using R ; Second Edition
#' @examples
#' ## We see that the number of repetitions n needs to be very large
#' ## in order to get even just two decimal places of accuracy.
#' f <- function(x) return(x^3 - 7 * x^2 + 1 )
#' hit_miss(f, 0, 1, -6, 2, 1000)
#' hit_miss(f, 0, 1, -6, 2, 10000)
#' hit_miss(f, 0, 1, -6, 2, 100000)
#' hit_miss(f, 0, 1, -6, 2, 1000000)


hit_miss <- function(ftn, a, b, f.min, f.max, n)
{
  Z.sum = 0
  for (i in 1:n)
  {
    X = runif(1, a, b)
    Y = runif(1, f.min, f.max)
    Z = (ftn(X) >= Y)
    Z.sum = Z.sum + Z
  }
  Total <- (b - a) * f.min + (Z.sum/n) * (b - a) * (f.max - f.min)
  return(Total)
}

#' Plot Monte-Carlo integration using the hit and miss method
#' @description  Partially vectorised version.
#' @param ftn a function of single variable.
#' @param a lower bound of the integration. Must be finite.
#' @param b upper bound of the integration. Must be finite and assume b > a.
#' @param c lower bound of the ftn over the range [a,b].
#' @param d upper bound of the ftn over the range [a,b].
#' @param n number of samples used in the estimation.
#' @export
#' @seealso \code{plot}
#' @references Owen Jones,Robert Maillardet,Andrew Robinson (2014).Introduction to Scientific Programming and Simulation Using R ; Second Edition
#' @examples
#' ## We have added a line to plot the successive approximations to the integral.
#' f <- function(x) return(x^3 - 7 * x^2 + 1 )
#' hmplot(f, 0, 1, -6, 2, 10000)
#' lines(c(0, 10000), c(-13/12, -13/12))


hmplot <- function(ftn, a, b, c, d, n)
{
  X <- runif(n, a, b)
  Y <- runif(n, c, d)
  Z <- (sapply(X, ftn) >= Y)
  Total <- (b - a)*c + (cumsum(Z)/(1:n)) * (b - a) * (d - c)
  plot(1:n, Total, type = "l",xlab = "Number of Points" , ylab = "Approximation of Integral")
  return(Total[n])
}

#' (Improved) Monte Carlo integration
#' @description Monte Carlo integral of ftn over [a, b] using a sample of size n.
#' @param a lower bound of the integration. Must be finite.
#' @param b upper bound of the integration. Must be finite and assume b > a.
#' @param n number of samples used in the estimation.
#' @details When we say a Monte Carlo technique is better than another we mean that using the same number of function calls, it has smaller variance.Because our estimates are based on random samples they are themselves random variables.
#' @export
#' @seealso \code{integral}
#' @references Owen Jones,Robert Maillardet,Andrew Robinson (2014).Introduction to Scientific Programming and Simulation Using R ; Second Edition
#' @examples
#' f <- function(x) return(x^3 + 1)
#' mc_integral(f, 0, 1, 10000)


mc_integral <- function(ftn, a, b, n)
{
  u <- runif(n, a, b)
  x <- sapply(u, ftn)
  return( mean(x) * (b-a) )
}


