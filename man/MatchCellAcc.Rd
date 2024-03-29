% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_cell.R
\name{MatchCellAcc}
\alias{MatchCellAcc}
\title{Match cell line name with cellosaurus accession}
\usage{
MatchCellAcc(names, file)
}
\arguments{
\item{names}{A vector of charactors contains names of interested cell lines}

\item{file}{The path to Cellosaurus XML file which was posted on
\url{https://web.expasy.org/cellosaurus/}}
}
\value{
A data frame contains:
\itemize{
  \item \strong{input_name} The input cell lines names
  \item \strong{cellosaurus_accession} The cellosaurus accession matched with
  cell line names.
  \item \strong{all_name} All names for the mached cell lines, including
  primary name and synonyms.
}
}
\description{
Match cell line name with cellosaurus accession
}
\author{
Jing Tang \email{jing.tang@helsinki.fi}
Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
}
