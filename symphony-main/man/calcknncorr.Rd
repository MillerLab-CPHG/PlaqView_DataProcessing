% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knncorr.R
\name{calcknncorr}
\alias{calcknncorr}
\title{Calculates the k-NN correlation, which measures how well the sorted ordering of k nearest reference
neighbors in a gold standard embedding correlate with the ordering for the same reference cells in
an alternative embedding (i.e. from reference mapping).
NOTE: it is very important for the order of reference cells (cols) in gold_ref matches that of alt_ref
(same for matching columns of gold_query and alt_query).}
\usage{
calcknncorr(gold_ref, alt_ref, gold_query, alt_query, k = 500)
}
\arguments{
\item{gold_ref}{Reference cells in gold standard embedding (PCs by cells)}

\item{alt_ref}{Reference cells in alternative embedding (PCs by cells)}

\item{gold_query}{Query cells in gold standard embedding (PCs by cells)}

\item{alt_query}{Query cells in alternative embedding (PCs by cells)}

\item{k}{Number of reference neighbors to use for kNN-correlation calculation}
}
\value{
Vector of k-NN correlations for query cells
}
\description{
Calculates the k-NN correlation, which measures how well the sorted ordering of k nearest reference
neighbors in a gold standard embedding correlate with the ordering for the same reference cells in
an alternative embedding (i.e. from reference mapping).
NOTE: it is very important for the order of reference cells (cols) in gold_ref matches that of alt_ref
(same for matching columns of gold_query and alt_query).
}
