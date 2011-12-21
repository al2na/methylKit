\docType{methods}
\name{read}
\alias{read}
\alias{read,character,character,character-method}
\alias{read,list,list,character-method}
\title{read file(s) to a methylrawList or methylraw object}
\usage{
  read(location,sample.id,assembly,pipeline="amp",header=T,
  context="CpG",treatment)
}
\arguments{
  \item{location}{file location(s), either a list of
  locations (each a character string) or one location
  string}

  \item{sample.id}{sample.id(s)}

  \item{assembly}{a string that defines the genome assembly
  such as hg18, mm9}

  \item{header}{if the input file has a header or not
  (default: TRUE)}

  \item{pipeline}{name of the alignment pipeline, currently
  only supports AMP (default: AMP)}

  \item{treatment}{a vector contatining 0 and 1 denoting
  which samples are control which samples are test}

  \item{context}{methylation context string, ex:
  CpG,CpH,CHH, etc. (default:CpG)}
}
\value{
  returns methylRaw or methylRawList
}
\description{
  read a list of locations or one location and create a
  methylrawList or methylraw object
}

