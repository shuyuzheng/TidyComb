% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_cell.R
\name{AnnotateCell}
\alias{AnnotateCell}
\title{Annotate cell lien from cell line name}
\usage{
AnnotateCell(cell_names, file)
}
\arguments{
\item{file}{The path to Cellosaurus XML file which was posted on
\url{https://web.expasy.org/cellosaurus/}}

\item{cell}{A vector of characters contains cell line names.}
}
\value{
A list with 2 data frames:
\itemize{
  \item \strong{cell_line} The cell_line table prepared for uploading.
  \item \strong{cell_id} The DrugComb IDs for new cell lines.
}
}
\description{
Annotate cell lien from cell line name
}
\author{
Jing Tang \email{jing.tang@helsinki.fi}
Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
}
