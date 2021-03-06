% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/yet_another_functions_up_to_date_24112018.R
\name{CreateDeltaAIPairwiseDF}
\alias{CreateDeltaAIPairwiseDF}
\title{Creates a table of parwise comparisons for techreps in given table.}
\usage{
CreateDeltaAIPairwiseDF(df, thrs = 2^c(0:12), thrsSide = "both",
  mlns = F, repcolmns = F, what = "noname")
}
\arguments{
\item{df}{A dataframe of genes/transcripts and parental counts for technical replicates in columns.}

\item{thrs}{An optional vector of thresholds (default = 2**c(0:12)).}

\item{thrsSide}{An optional parameter of threshold side: 'both' or 'low' (default = 'both').}

\item{mlns}{An optionsl binary parameter: if the file contains data of millionr read sampling, FALSE or TRUE (default = F).}

\item{repcolmns}{An optional parameter for a range op replicates for consideration (default = all replicates in df).}

\item{what}{A name, is needed if not mlns and no names in list (default = "noname").}
}
\value{
A table of parwise comparisons for techreps in given table.
}
\description{
Creates a table of parwise comparisons for techreps in given table.
}
\examples{

}
