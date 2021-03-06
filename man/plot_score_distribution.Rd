% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_plots.R
\name{plot_score_distribution}
\alias{plot_score_distribution}
\title{Drawing a score distribution plot}
\usage{
plot_score_distribution(traceQC_input)
}
\arguments{
\item{traceQC_input}{A TraceQC object}
}
\value{
A ggplot2 object that shows alignment score distribution.
}
\description{
Drawing a score distribution plot
}
\examples{
data(example_obj)
plot_score_distribution(example_obj)

}
