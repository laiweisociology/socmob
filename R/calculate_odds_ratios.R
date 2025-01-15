#' Calculate Odds Ratios for a Contingency Table
#'
#' This function computes all possible odds ratios from a given contingency table.
#' For 2x2 tables, it returns a single odds ratio.
#' For larger tables, it calculates the odds ratio for each possible 2x2 subtable.
#'
#' @param tbl A contingency table (matrix or table) with at least two rows and two columns.
#'
#' @return
#' - For 2x2 tables: A single numeric value representing the odds ratio.
#' - For larger tables: A named numeric vector containing the odds ratios for each 2x2 subtable.
#'   Odds ratios involving zero cells are returned as `NA`.
#'
#' @examples
#' # Example with a 2x2 table
#' tbl_2x2 <- matrix(c(10, 20, 30, 40), nrow = 2, byrow = TRUE)
#' calculate_odds_ratios(tbl_2x2)
#'
#' # Example with a 3x3 table
#' tbl_3x3 <- matrix(c(10, 20, 30,
#'                    20, 30, 40,
#'                    30, 40, 50), nrow = 3, byrow = TRUE)
#' calculate_odds_ratios(tbl_3x3)
#'
#' @export
calculate_odds_ratios <- function(tbl) {
  # Validate input
  if (!is.matrix(tbl) && !is.table(tbl)) {
    stop("Input must be a matrix or table.")
  }

  nr <- nrow(tbl)
  nc <- ncol(tbl)

  if (nr < 2 || nc < 2) {
    stop("Table must have at least two rows and two columns.")
  }

  # Function to compute odds ratio for a 2x2 subtable
  compute_or <- function(sub_tbl) {
    a <- sub_tbl[1,1]
    b <- sub_tbl[1,2]
    c <- sub_tbl[2,1]
    d <- sub_tbl[2,2]
    if (b * c == 0) {
      return(NA)  # Undefined odds ratio due to zero cell(s)
    }
    return((a * d) / (b * c))
  }

  # If table is 2x2, compute and return the single odds ratio
  if (nr == 2 && nc == 2) {
    return(compute_or(tbl))
  }

  # For larger tables, compute odds ratios for all 2x2 subtables
  row_combs <- combn(nr, 2)
  col_combs <- combn(nc, 2)

  odds_ratios <- c()

  for (i in 1:ncol(row_combs)) {
    for (j in 1:ncol(col_combs)) {
      rows <- row_combs[,i]
      cols <- col_combs[,j]
      sub_tbl <- tbl[rows, cols]
      or <- compute_or(sub_tbl)
      name <- paste("Rows", rows[1], "&", rows[2],
                    "Cols", cols[1], "&", cols[2], sep = "_")
      odds_ratios[name] <- or
    }
  }

return(odds_ratios)
}
