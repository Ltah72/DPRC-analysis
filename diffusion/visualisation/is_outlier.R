#calculate outliers and visualise on graph

#source: https://stackoverflow.com/questions/33524669/labeling-outliers-of-boxplots-in-r

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}