\docType{methods}
\name{read}
\alias{read}
\alias{read,-methods}
\alias{read,character,character,character,ANY-method}
\alias{read,list,list,character,ANY-method}
\alias{read,methylRawList-method}
\title{read file(s) to a methylrawList or methylraw object}
\arguments{
  \item{file}{location(s)}

  \item{sample.id}{sample.id(s)}

  \item{assembly}{a string that defines the genome assembly
  such as hg18, mm9}

  \item{header}{if the input file has a header or not
  (default: TRUE)}

  \item{pipeline}{name of the alignment pipeline, currently
  only supports AMP (default: AMP)}

  \item{treatment}{a vector contatining 0 and 1 denoting
  which samples are control which samples are test}
}
\description{
  read a list of locations or one location and create a
  methylrawList or methylraw object
}

