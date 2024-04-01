# ----------------------- Helper Functions to Implement ------------------------

#' Evaluate whether the argument is less than 2
#'
#' Returns TRUE if the numeric argument x is a prime number, otherwise returns
#' FALSE
#'
#' @param x (numeric): the numeric value(s) to test
#'
#' @return logical value or vector indicating whether the numeric argument is less than 2
#' @export
#'
#' @examples
#' less_than_zero(-1)
#' [1] TRUE
#' less_than_zero(10)
#' [1] FALSE
#' less_than_zero(c(-1,0,1,2,3,4))
#' [1] TRUE FALSE FALSE FALSE FALSE FALSE
less_than_zero <- function(x) {
  x < 0
}

#' Evaluate whether the argument is between two numbers
#'
#' Returns TRUE if the numeric argument x is contained within the closed interval
#' (a, b), otherwise return FALSE.
#'
#' @param x (numeric): the numeric value(s) to test
#' @param a (number): the lower bound
#' @param b (number): the upper bound
#'
#' @return logical value of same type as input (e.g. scalar, vector, or matrix)
#' @export
#'
#' @examples
#' is_between(3,1,5)
#' [1] TRUE
#' is_between(c(1,9,5,2), 1, 5)
#' [1] FALSE FALSE FALSE TRUE
#' is_between(matrix(1:9, nrow=3, byrow=TRUE), 1, 5)
#'       [,1]  [,2]  [,3]
#' [1,] FALSE  TRUE  TRUE
#' [2,]  TRUE FALSE FALSE
#' [3,] FALSE FALSE FALSE
is_between <- function(x, a, b) {
  x > a & x < b
}

#' Return the values of the input vector that are not NA
#'
#' Returns the values of the input vector `x` that are not NA
#'
#' @param x (numeric): numeric vector to remove NA from 
#'
#' @return numeric vector `x` with all `NA` removed
#' @export
#'
#' @examples
#' x <- c(1,2,NA,3)
#' rm_na(x)
#' [1] 1 2 3
rm_na <- function(x) {
  x <- x[!is.na(x)]
  return(x)
}

#' Calculate the minimum of each row of a matrix
#'
#' Given the matrix x with n rows and m columns, return a numeric vector of
#' length n that contains the minimum value of each row of x
#'
#' @param x (numeric matrix): matrix to compute minumum along rows
#'
#' @return (numeric vector) vector containing minimum row values
#' @export
#'
#' @examples
#' m <- matrix(1:9, nrow=3, byrow=T)
#' row_mins(m)
#' [1] 1 4 7
#' 
row_medians <- function(m) {
  apply(m, MARGIN = 1, median, na.rm=TRUE)
}

#' Evaluate each row of a matrix with a provided function
#'
#' Given the matrix `x` with n rows and m columns, return a numeric vector of
#' length n that contains the returned values from the evaluation of the
#' function `fn` on each row. `fn` should be a function that accepts a vector as
#' input and returns a scalar value
#'
#' @param x (numeric matrix): matrix to evaluate the function along rows
#' @param fn (function) function that accepts a vector as input and returns a scalar
#' @param na.rm (logical) OPTIONAL: a logical evaluating to `TRUE` or `FALSE`
#'   indicating whether `NA` values should be stripped before computing summary
#'
#' @return (numeric vector) vector of the evaluation of `fn` on rows of `x`
#' @export
#'
#' @examples
#' m <- matrix(1:9, nrow=3, byrow=T)
#' summarize_rows(m, min)
#' [1] 1 4 7
#' summarize_rows(m, mean)
#' [1] 2 5 8
summarize_rows <- function(x, fn, na.rm=FALSE) {
  apply(x, MARGIN = 1, fn)
}

#' Summarize matrix rows into data frame
#'
#' Summarizes the rows of matrix `x`. Returns a data frame with the following
#' columns in order and the corresponding value for each row:
#' 
#'   * mean - arithmetic mean
#'   * stdev - standard deviation
#'   * median - median
#'   * min - minimum value
#'   * max - maximum value
#'   * num_lt_0 - the number of values less than 0
#'   * num_btw_1_and_5 - the number of values between 1 and 5
#'   * num_na - the number of missing (NA) values
#'
#' @param x (numeric matrix): matrix to evaluate the function along rows
#' @param na.rm (logical) OPTIONAL: a logical evaluating to `TRUE` or `FALSE`
#'   indicating whether `NA` values should be stripped before computing summary
#'
#' @return (data frame) data frame containing summarized values for each row of `x`
#' @export
#'
#' @examples
#' m <- matrix(1:9, nrow=3, byrow=T)
#' summarize_matrix(m)
#'   mean stdev median min max num_lt_0 num_btw_1_and_5 num_na
#' 1    2     1      2   1   3        0               3      0
#' 2    5     1      5   4   6        0               1      0
#' 3    8     1      8   7   9        0               0      0
#'
m <- matrix(rnorm(1000), nrow=4, byrow=T)
summarize_matrix(m)
#'          mean    stdev      median       min      max num_lt_0 num_btw_1_and_5 num_na
#' 1 -0.02220010 1.006901  0.02804177 -3.485147 3.221089      120              61      0
#' 2 -0.01574033 1.026951 -0.04725656 -2.967057 2.571608      112              70      0
#' 3 -0.09040182 1.027559 -0.02774705 -3.026888 2.353087      130              54      0
#' 4  0.09518138 1.030461  0.11294781 -3.409049 2.544992       90              72      0
summarize_matrix <- function(x, na.rm = FALSE) {
  means <- apply(x, MARGIN = 1, mean)
  st_dev <- apply(x, MARGIN = 1, sd)
  medians <- apply(x, MARGIN = 1, median)
  mins <- apply(x, MARGIN = 1, min)
  maxes <- apply(x, MARGIN = 1, max)
  num_less_than_zero <- sum(summarize_rows(x, less_than_zero))
  num_between_1_and_5 <- apply(is_between(x, 1, 5), 1, sum)
  num_na <- apply(is.na(x), MARGIN = 1, sum)
  
  result_dataframe <- data.frame(mean = means, stdev = st_dev, median = medians,
                                 min = mins, max = maxes, num_lt_0 = num_less_than_zero,
                                 num_btw_1_and_5 = num_between_1_and_5, num_na = num_na)
  return(result_dataframe)
}

# ------------ Helper Functions Used By Assignment, You May Ignore ------------
sample_normal <- function(n, mean=0, sd=1) {
  set.seed(1337)
  samples <- rnorm(n, mean=mean, sd=sd)
  return(samples)
}

sample_normal_w_missing <- function(n, mean=0, sd=1, missing_frac=0.1) {
  set.seed(1337)
  samples <- rnorm(n, mean=mean, sd=sd)
  missing <- rbinom(length(samples), 1, missing_frac)==1
  samples[missing] <- NA
  return(samples)
}

simulate_gene_expression <- function(num_samples, num_genes) {
  set.seed(1337)
  gene_exp <- matrix(
    rnbinom(num_samples*num_genes, rlnorm(num_genes,meanlog = 3), prob=runif(num_genes)),
    nrow=num_genes
  )
  return(gene_exp)
}

simulate_gene_expression_w_missing <- function(num_samples, num_genes, missing_frac=0.1) {
  gene_exp <- simulate_gene_expression(num_samples, num_genes)
  missing <- matrix(
    rbinom(num_samples*num_genes, 1, missing_frac)==1,
    nrow=num_genes
  )
  gene_exp[missing] <- NA
  return(gene_exp)
}