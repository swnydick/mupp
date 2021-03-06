% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_mupp.R
\name{simulate_mupp_resp}
\alias{simulate_mupp_resp}
\title{Simulate MUPP Responses}
\usage{
simulate_mupp_resp(persons, items)
}
\arguments{
\item{persons}{persons data.frame with column names
\verb{[person, dim, theta]}}

\item{items}{items data.frame with column names
\verb{[item, statement, dim, alpha, delta, tau]}}
}
\value{
a data.frame of \verb{[person x item x response pattern]}
}
\description{
Generate responses that can be used for the MUPP model.
}
\details{
The persons and items df needs to look identical to that coming from
\code{\link{simulate_mupp_params}} or else this function will not work.
}
\author{
Steven Nydick, \email{steven.nydick@kornferry.com}
}
