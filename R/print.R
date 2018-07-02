#' @export
#' @import methods
print.bayesdfa <- function(x, ...) {
  base::print(x$monitor, digits = 2)
}

# test
