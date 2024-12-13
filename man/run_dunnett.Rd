% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dunnett_test.R
\name{run_dunnett}
\alias{run_dunnett}
\title{Performs Dunnett's T3 test on each row of a data.frame}
\usage{
run_dunnett(my_data, var1)
}
\arguments{
\item{my_data}{a numeric data.frame}

\item{var1}{a list which maps columns of my_data to experimental factors}
}
\value{
a data.frame of p values
}
\description{
Performs Dunnett's T3 test on each row of a data.frame
}
