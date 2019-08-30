\name{LPIexpand.input}
\title{LPIexpand.input}
\alias{LPIexpand.input}
\description{Expand vectors or lists to the LPI internal format with
  names "RX1", "RX2", "TX1", and "TX2"}
\usage{LPIexpand.input( parvec)}

\arguments{
  \item{ parvec }{ A vector or list. See details.}
  }

  \details{
    The input may be in several formats:

    1. If parvec is readily in the internal format it is returned as
    such.

    2. If parvec does not have named elements, it is repeated /
    truncated to length 4 and the elements are named in order
    "RX1", "RX2", "TX1", "TX2".

    3. If parvec does not contain all the internally used names, but it
    contains also unnamed elements, the unnamed ones are repeated and
    named to fill the output.

    4. Names 'RX' and 'TX' are expanded in an obvious way.

    5. Names 'TR1' and 'TR2' are expanded in an obvious way.
  }
  \value{
    A named vector or list of the internally used format
  }

  

  \examples{

# Expansion of a single unnamed value
LPIexpand.input( 4 )

# A combination of named and unnamed elements is accepted
LPIexpand.input( c( RX1=2, 4) )

# The RX and TX entries
LPIexpand.input(c(RX=2,TX=3))

# The TR1 and TR2 entries
LPIexpand.input(c(TR1=2,TR2=3))

}



  \author{Ilkka Virtanen (University of Oulu, Finland) \cr
  \email{ilkka.i.virtanen@oulu.fi}}
