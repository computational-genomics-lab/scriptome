% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tr2g.R
\name{tr2g_TxDb}
\alias{tr2g_TxDb}
\title{Get transcript and gene info from TxDb objects}
\usage{
tr2g_TxDb(txdb)
}
\arguments{
\item{txdb}{A \code{\link{TxDb}} object with gene annotation.}
}
\value{
A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
for transcript ID, and \code{tx_id} for internal transcript IDs used to avoid
duplicate transcript names. For TxDb packages from Bioconductor, gene ID is
Entrez ID, while transcript IDs are Ensembl IDs with version numbers for
\code{TxDb.Hsapiens.UCSC.hg38.knownGene}. In some cases, the transcript ID
have duplicates, and this is resolved by adding numbers to make the IDs
unique.

A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
for transcript ID, and \code{gene_name} for gene names. If \code{other_attrs}
has been specified, then those will also be columns in the data frame returned.
}
\description{
The genome and gene annotations of some species can be conveniently obtained
from Bioconductor packages. This is more convenient than downloading GTF
files from Ensembl and reading it into R. In these packages, the gene
annotation is stored in a \code{\link{TxDb}} object, which has standardized
names for gene IDs, transcript IDs, exon IDs, and so on, which are stored in
the metadata fields in GTF and GFF3 files, which are not standardized.
This function extracts transcript and corresponding gene information from
gene annotation stored in a \code{\link{TxDb}} object.
}
\examples{
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
tr2g_TxDb(TxDb.Hsapiens.UCSC.hg38.knownGene)
}
\seealso{
Other functions to retrieve transcript and gene info: \code{\link{sort_tr2g}},
  \code{\link{tr2g_EnsDb}}, \code{\link{tr2g_ensembl}},
  \code{\link{tr2g_fasta}}, \code{\link{tr2g_gff3}},
  \code{\link{tr2g_gtf}}, \code{\link{transcript2gene}}

Other functions to retrieve transcript and gene info: \code{\link{sort_tr2g}},
  \code{\link{tr2g_EnsDb}}, \code{\link{tr2g_ensembl}},
  \code{\link{tr2g_fasta}}, \code{\link{tr2g_gff3}},
  \code{\link{tr2g_gtf}}, \code{\link{transcript2gene}}
}
\concept{functions to retrieve transcript and gene info}
