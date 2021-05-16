




#dependancy
require(stringr)

#' Get tair identifier
#' 
#' Extract arabidopsis (TAIR) gene identifiers from a text. These identifiers conform to the following structure:
#' AT followed by the chromosome number the gene is located on; 
#' or alternatively a M for mitochondria or C for chloroplast if the gene is located on one of the plastid genomes.
#' 
#' @param string The text containing the gene identifiers.
#' @return A list of Tair identifiers found in string.
#' @examples 
#' ids<-gettairid(string=c("AT1G19030 AT3G61800,AT2G02200convallisAT5G23670,ATMG00730loreAT3G28291,AT4G07460,AT3G29ATMG00730"))
#' 
#' 
#' 
gettairid<-function(string){
  require(stringr)
  ls_geneid<-str_extract_all(string = string,
                             pattern =  
                               "AT(M|C|[:digit:]{1})(G{1})([:digit:]{5})")
  
  return(ls_geneid)
}