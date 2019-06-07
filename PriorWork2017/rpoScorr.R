
# the correlate function from the corrr package
# I didn't want the diagonal to be NA

correlate2 <- function (x, y = NULL, use = "pairwise.complete.obs", method = "pearson") 
{
  x <- stats::cor(x = x, y = y, use = use, method = method)
  x <- dplyr::as_data_frame(x)
  x <- corrr::first_col(x, names(x))
  class(x) <- c("cor_df", class(x))
  x
}

