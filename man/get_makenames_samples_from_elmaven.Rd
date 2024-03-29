% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_makenames_samples_from_elmaven.R
\name{get_makenames_samples_from_elmaven}
\alias{get_makenames_samples_from_elmaven}
\title{get_makenames_samples_from_elmaven}
\usage{
get_makenames_samples_from_elmaven(maven_output_df = NULL)
}
\arguments{
\item{maven_output_df}{The peak detailed format of elmaven}
}
\value{
The named vector of samplenames using make.names
}
\description{
It returns a named vector of samplenames where make.names are assigned to samples.
}
\examples{
get_makenames_samples_from_elmaven(maven_output_df)
}
