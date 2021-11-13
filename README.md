# gaRdenbox

A collection of miscellanious R functions usefull to me. Currently included:

1. gettairid, a function to extract the first Arabidopsis gene ID from a character string (following the TAIR nomenclature). The function also allowes for the selection of specific chromosomes or plastids. Since this function does not return a list it is more suitable for dataframe operations.

2. gettairid_all, a function to extract Arabidopsis gene ID's from a character string (following the TAIR nomenclature). The function also allowes for the selection of specific chromosomes or plastids. Since this function returns a list of all geneid's in the query string it is more suitable for text extraction.


You are welcome to use them!

To install the package please run:

>install.packages("devtools") 

>devtools::install_github("KTMD-plant/gaRdenbox")
