% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_group_summary_from_peak_detailed.R
\name{make_group_summary_from_peak_detailed}
\alias{make_group_summary_from_peak_detailed}
\title{make_group_summary_from_peak_detailed}
\usage{
make_group_summary_from_peak_detailed(
  maven_output_df = NULL,
  quant_type = "peakAreaCorrected"
)
}
\arguments{
\item{maven_output_df}{Provide peak detailed format of elmaven}

\item{quant_type}{The type of qunatification to be used which is from peak detailed format columns}
}
\value{
A dataframe of group summary format
}
\description{
It converts the peak detailed to group summary format of elmaven
}
