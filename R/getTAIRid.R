




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


#function to get all tair id's regardless of chromosme or plastid
f_allid<-function(querytext){
  ls_geneid<-stringr::str_extract_all(string = querytext,
                                      pattern =  
                                        "AT(M|C|[:digit:]{1})(G{1})([:digit:]{5})")
  
  return(ls_geneid)
  
}

#function to only get nucleaic genes
f_nuclid<-function(querytext){
  ls_geneid<-stringr::str_extract_all(string = querytext,
                                      pattern =  
                                        "AT([:digit:]{1})(G{1})([:digit:]{5})")
  
  return(ls_geneid)
  
}

#function to only get plastid genes
f_plastid<-function(querytext){
  ls_geneid<-stringr::str_extract_all(string = querytext,
                                      pattern =  
                                        "AT(M|C{1})(G{1})([:digit:]{5})")
  
  return(ls_geneid)
  
}

# #function to validate exclusion input
# f_validate_excl<-function(exclusionstatement)

#function to get specified chromosomes or M or C

f_spec<-function(querytext, specify, validation){
  
  #get specs as reg expr
  specs<-paste(specify, collapse = "|")
  
  ls_geneid<-stringr::str_extract_all(string = querytext,
                                      pattern =
                                        paste0("AT([", specs, "]{1})(G{1})([:digit:]{5})")
                                      )
  return(ls_geneid)
}

#function checking validity of specified vector

f_validate_spec<-function(specify){
  
  
  ##check if specify is a vector
  if (!is.vector(specify)) {
    stop("Error: specify must be a character vector")
    
  }
  #check length of vector >0
  if (!length(specify)>=1) {
    stop("Error: specify must be a vector of length 1 or more")
    
  }
  #check if specify contains characters
  if (!is.character(specify)) {
    stop("Error: specify must be a charcter vector. Please put chromosome numbers in quotes.")
    
  }
  #allowed variables in the specify vector:
  allowed<-c("1","2","3","4","5","M","C")
  
  #contains forbidden values
  
  forbidden<-!specify%in%allowed
  if (TRUE%in%forbidden) {
    stop("Error: please only specify chromosome numbers and/or the letters \"M\" and \"C\" ")
    
  }
  
  
  return(TRUE)
}


#MAIN
gettairid<-function(querytext, exclude_plastid, exclude_nuclear, specify) {                                               
  #if all optional variables are left empty call allid
  if(
    missing(exclude_plastid) & 
    missing(exclude_nuclear) & 
    missing(specify)){
    ls_geneid<-f_allid(querytext = querytext)
  }
  
  return(ls_geneid)
  
  #if exclude plastid is called check for FALSE >pass on, TRUE >call f_nuclid, or ELSE >raise error
  if(
    !missing(exclude_plastid) & 
    missing(exclude_nuclear) & 
    missing(specify)){
    if (exclude_plastid==FALSE) {
      explas=FALSE
    }
    if (exclude_plastid==TRUE) {
      explas=TRUE
      #call f_nuclid
      ls_geneid<-f_nuclid(querytext = querytext)
      
      return(ls_geneid)
    }
    else{
      stop("Error: exclude_plastid is an optional variable. If used it needs to be TRUE or FALSE.")
    }
  }
  
  
  #if exclude nuclear is called check for FALSE >pass on, TRUE >call f_plastid, or ESLE >raise error
  if(
    missing(exclude_plastid) & 
    !missing(exclude_nuclear) & 
    missing(specify)){
    if (exclude_nuclear==FALSE) {
      exnucl=FALSE
    }
    if (exclude_nuclear==TRUE) {
      exnucl=TRUE
      #call f_plastid function
      ls_geneid<-f_plastid(querytext = querytext)
      
      return(ls_geneid)
    }
    else{
      stop("Error: exclude_nuclear is an optional variable. If used it needs to be TRUE or FALSE.")
    }
  }
  
  #if specifiy is called check if it is a vector
  if(
    missing(exclude_plastid) & 
    missing(exclude_nuclear) & 
    !missing(specify)){
    
    #validate specified
    pass<-f_validate_spec(specify = specify)
    
    #call f_spec
  ls_geneid<-f_spec(querytext = querytext,specify = specify, validation = pass)
  
  return(ls_geneid)
  }
  
  
  ##edge cases
  #if exclude nuclear is true and specify is spedified >raise error
 
  
  
  #if explas & exnucl ==FALSE raise ERROR, wtf do you want? 


}
