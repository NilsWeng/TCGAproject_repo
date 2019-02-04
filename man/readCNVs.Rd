\name{readCNVs}
\alias{readCNVs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Given a filename of a TCGA cnv file opens and convert it to GRanges object
}
\description{
Given a filename of a TCGA cnv file opens and convert it to GRanges object
}
\usage{
readCNVs(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
%%     ~~Describe \code{filename} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (filename) 
{
    suppressPackageStartupMessages(library("GenomicRanges"))
    cnv_DF <- read.table(filename, header = TRUE, sep = "\t", 
        stringsAsFactors = FALSE)
    cnv_GRanges <- GRanges(seqnames = cnv_DF$Chromosome, ranges = IRanges(cnv_DF$Start, 
        cnv_DF$End), strand = "*", Num_Probes = cnv_DF$Num_Probes, 
        Segment_Mean = cnv_DF$Segment_Mean,
        GDC_Aliquot = cnv_DF$GDC_Aliquot)
        
    return(cnv_GRanges)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
