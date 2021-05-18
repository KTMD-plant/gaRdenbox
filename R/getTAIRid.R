




#dependancy

#' Get TAIR identifier
#' 
#' Extract arabidopsis (TAIR) gene identifiers from a text. These identifiers conform to the following structure:
#' 
#' "AT" 
#' followed by the chromosome number the gene is located on; 
#' 
#' or alternatively a M for mitochondria or C for chloroplast if the gene is located on one of the plastid genomes.
#' 
#' "G" followed by 
#' five-digit code, numbered from top/north to bottom/south of chromosome
#' 
#' Example of gene id: 
#' nuclear: AT3G28291
#' mitochondrial: ATMG00730
#' 
#' Does not extract: transposable elements (these do not follow this naming convention)
#' For more info on TAIR ID nomenclature see https://www.arabidopsis.org/portals/nomenclature/guidelines.jsp
#' 
#' 

#' 
#' @param querytext The text (character string) containing the gene identifiers. Make sure the character string does not contain quotation marks within.
#' @return A list of TAIR identifiers found in the querytext.
#' @examples 
#' ids<-gettairid(querytext="AT1G19030 AT3G61800,AT2G02200convallisAT5G23670,ATMG00730loreAT3G28291,AT4G07460,AT3G29ATMG00730"
#' 
#' 
#' 
#' 
gettairid<-function(querytext){
  ls_geneid<-stringr::str_extract_all(string = querytext,
                             pattern =  
                               "AT(M|C|[:digit:]{1})(G{1})([:digit:]{5})")
  
  return(ls_geneid)
}
