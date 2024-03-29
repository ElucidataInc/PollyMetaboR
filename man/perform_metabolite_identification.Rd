% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perform_metabolite_identification.R
\name{perform_metabolite_identification}
\alias{perform_metabolite_identification}
\title{perform_metabolite_identification}
\usage{
perform_metabolite_identification(
  mz_data = NULL,
  comp_data = NULL,
  mz_colname = "basemass",
  mz_tolerence_unit = "ppm",
  mz_tolerence = 20,
  rt_tolerence = NULL,
  rt_colname = "rt",
  wrap_comp = FALSE,
  numcores = 2
)
}
\arguments{
\item{mz_data}{A dataframe or numeric vector of mz values}

\item{comp_data}{The compound database used for identification}

\item{mz_colname}{The mz column name present in mz_data dataframe}

\item{mz_tolerence_unit}{The mz tolerance unit (ppm or Da)}

\item{mz_tolerence}{Value of mz tolerence}

\item{rt_tolerence}{Value of rt tolerence in minutes}

\item{rt_colname}{The rt column name present in mz_data dataframe}

\item{wrap_comp}{Wrap compounds in same row detected for same feature (mz and rt)}

\item{numcores}{Number of cores used for processing}
}
\value{
A dataframe with identified metabolites
}
\description{
It restructures the camera output to group features
}
\examples{
perform_metabolite_identification(mz_data, comp_data, mz_colname = 'basemass',
                                  mz_tolerence_unit = "ppm", mz_tolerence = 20)
}
