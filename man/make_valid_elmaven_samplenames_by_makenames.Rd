% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_valid_elmaven_samplenames_by_makenames.R
\name{make_valid_elmaven_samplenames_by_makenames}
\alias{make_valid_elmaven_samplenames_by_makenames}
\title{make_valid_elmaven_samplenames_by_makenames}
\usage{
make_valid_elmaven_samplenames_by_makenames(
  maven_output_df = NULL,
  metadata_df = NULL
)
}
\arguments{
\item{maven_output_df}{The peak detailed format of elmaven}

\item{metadata_df}{The sample to cohort mapping file where samples should be present in first column}
}
\value{
The list of two dataframes with replaced samples names by make.names
}
\description{
It makes valid samplenames in elmaven file using make names.
}
\examples{
make_valid_elmaven_samplenames_by_makenames(maven_output_df, metadata_df)
}
