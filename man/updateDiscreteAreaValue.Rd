\name{updateDiscreteAreaValue}
\alias{updateDiscreteAreaValue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Update value of discreteArea object.
}
\description{
Update value(s) of one or some areas of a discreteArea object.
}
\usage{
updateDiscreteAreaValue(x, newval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a discreteArea object: list(areaValues, discretePoints), where
    areaValues: data.frame(areaId,centx,centy,value)
    discretePoints: data.frame(areaId,ptx,pty,weight)
}
  \item{newval}{
  new values: a dataframe(areaId, value).
  }
}
\value{
a new discreteArea.
}
