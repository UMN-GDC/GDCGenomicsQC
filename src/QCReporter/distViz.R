library(tidyverse)


#' @title Histogram plotter'
#' @description This function plots a histogram of a column in a Dataframe
#' @param df A DataFrame
#' @param colName A column name in the DataFrame
#' @param cutoff A cutoff value for the histogram
#' @param title A title for the histogram
#' @param xaxislabel A label for the x-axis
#' @return A histogram plot
#' @examples
#' histPlotter(df, MAF, 0.05, "MAF", "MAF")
histPlotter <- function(df, colName, cutoff, title, xaxislabel){
  inflateds <- df %>% 
    summarize(zeros = sum({{ colName }} == 0), ones = sum({{ colName }} == 1))
  
  df %>% 
    filter({{ colName }} != 0 & {{ colName }} != 1) %>% 
    ggplot(aes(x = {{ colName }})) +
    geom_histogram() +
    geom_vline(xintercept = cutoff, color = "red") +
    annotate("text", x = -Inf, y = Inf, label = paste("zeros:", inflateds$zeros), vjust = 1, hjust = 0) +
    annotate("text", x = Inf, y = Inf, label = paste("ones:", inflateds$ones), vjust = 1, hjust = 1) +
    ggtitle(title) +
    xlab(xaxislabel) +
    theme_bw()
}

