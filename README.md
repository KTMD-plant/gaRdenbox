# gaRdenbox

A collection of miscellanious R functions usefull to me. Currently included:

1. gettairid, a function to extract the first Arabidopsis gene ID from a character string (following the TAIR nomenclature). The function also allowes for the selection of specific chromosomes or plastids. Since this function does not return a list it is more suitable for dataframe operations.
2. gettairid_all, a function to extract Arabidopsis gene ID's from a character string (following the TAIR nomenclature). The function also allowes for the selection of specific chromosomes or plastids. Since this function returns a list of all geneid's in the query string it is more suitable for text extraction.
3. goseq_nm, a function performing Gene Ontology enrichment for non-model species.
4. get_peaks a function to extract peaks from Malley et al 2016 (doi: 10.1016/j.cell.2016.04.038) of  a given (promotor) section of interest.
5. get_tf_from_peaks to get unique geneids from transcription factors binding in given promotor section.


You are welcome to use them!

To install the package please run:

>install.packages("devtools") 

>devtools::install_github("KTMD-plant/gaRdenbox")
