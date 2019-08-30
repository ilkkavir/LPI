\name{stripACF}
\title{stripACF}
\alias{stripACF}
\description{Strip unwanted ranges and lags from ACFs}
\usage{stripACF(ACFlist, rgates, lags, fullCovar = FALSE)}

\arguments{
  \item{ ACFlist }{ An ACF list read from an LPI output file or returned
    by LPIrun. }
  \item{ rgates }{ Indices of range gates to preserve }
  \item{ lags }{ Indices of lags to preserve }
  \item{ fullCovar }{ Logical, TRUE if the ACF list contains full
  covariance matrices. }

  }

  \value{
    An ACF list with only the selected ranges and lags maintained.
  }


  \examples{
\dontrun{

  # Load an LPI result file
  load( "1345541810000LP.Rdata" )

  # There will now be an ACF list on the workspace.
  # Check if it contains a full covariance matrix
  fullCov <- ifelse( is.null(ACF$covariance) , FALSE , TRUE )
  # Select range gates 30 - 40 and lags 1, 4, and 5.
  ACFstrip <- stripACF( ACF , seq(30,40) , c(1,4,5) , fullCov )
}
}



  \author{Ilkka Virtanen (University of Oulu, Finland) \cr
  \email{ilkka.i.virtanen@oulu.fi}}
