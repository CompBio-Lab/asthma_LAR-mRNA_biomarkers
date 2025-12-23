#' customeTheme function for ggplot
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param sizeStripFont font of size of facet labels
#' @param xAngle angle of x-axis labels
#' @param hjust horizontal justification 0-left, 0.5-center, 1-right
#' @param vjust vertical justification 0-low, 0.5-middle, 1-high
#' @param xSize font size of x-axis label
#' @param ySize font size of y-axis label
#' @param xAxisSize font size of x-axis label title
#' @param yAxisSize fotn size of y-axis label title
#' @export
customTheme = function(sizeStripFont, xAngle, hjust, vjust, xSize,
                       ySize, xAxisSize, yAxisSize) {
  theme(strip.background = element_rect(colour = "black", fill = "white",
                                        size = 1), strip.text.x = element_text(size = sizeStripFont),
        strip.text.y = element_text(size = sizeStripFont), axis.text.x = element_text(angle = xAngle,
                                                                                      hjust = hjust, vjust = vjust, size = xSize, color = "black"),
        axis.text.y = element_text(size = ySize, color = "black"),
        axis.title.x = element_text(size = xAxisSize, color = "black"),
        axis.title.y = element_text(size = yAxisSize, color = "black"),
        panel.background = element_rect(fill = "white", color = "black"))
}