



#FUNCTIONS USED IN GETTAIRID



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

f_logictest<-function(logic_var, name_var){

  if (!is.logical(logic_var)) {
    stop(
      paste("Logic error: ", name_var, " exclude_nuclear must be TRUE, FALSE or left empty.")
    )
  }

}

#MAIN

#dependancy

#' @title  Get TAIR identifier
#' 
#' @description 
#' Extract arabidopsis (TAIR) gene identifiers from a text. 
#' 
#' @details 
#' These identifiers conform to the following structure:
#' 
#' "AT" followed by the 
#' 
#' Chromosome number, a "M" for mitochondria or a "C" for chloroplast.
#' 
#' "G" followed by a
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
#' @param exclude_plastid OPTIONAL - When TRUE excludes identifiers of plastid genes from output (mitochondrial or chloroplastic).
#' @param exclude_nuclear OPTIONAL - When TRUE excludes identifiers of nuclear genes from output.
#' @param specify OPTIONAL - Allows you to specify chromosome numbers or "M" or "C" to be included in output.
#' Takes character vector of chromosome numbers and/or "M" and "C" for mitochondrial or chloroplastic gene identifiers respectively.
#' remember to put chromosome numbers in quotation marks.
#' 
#' 
#' When OPTIONAL variables are left empty all gene identifiers are extracted.
#' 
#' @return A list of TAIR identifiers found in the querytext.
#' @examples 
#' q<-"AT1G19030 AT3G61800,AT2G02200convallisAT5G23670,ATMG00730loreAT3G28291,AT4G07460,AT3G29ATMG00730"
#' 
#' ids<-gettairid(querytext=q)
#' ids:
#' "AT1G19030" "AT3G61800" "AT2G02200" "AT5G23670" "ATMG00730" "AT3G28291" "AT4G07460" "ATMG00730"
#' 
#' ids<-gettairid(querytext = q, exclude_nuclear = TRUE)
#' ids: "ATMG00730" "ATMG00730"
#' 
#' ids<-gettairid(querytext = q, specify = c("3", "5"))
#' ids: "AT3G61800" "AT5G23670" "AT3G28291"
#' 
#' 
#' 
#' @author Kilian Duijts
#' 
#' 
#' 
#' @export
gettairid<-function(querytext, exclude_plastid, exclude_nuclear, specify) {                                               
  #if all optional variables are left empty call allid
  if(
    missing(exclude_plastid) & 
    missing(exclude_nuclear) & 
    missing(specify)){
    ls_geneid<-f_allid(querytext = querytext)
  }
  

  #if exclude plastid is called check for FALSE >pass on, TRUE >call f_nuclid, or ELSE >raise error
  if(
    !missing(exclude_plastid) & 
    missing(exclude_nuclear) & 
    missing(specify)){
    
    #call logic test
    f_logictest(logic_var = exclude_plastid, name_var = "exclude_plastid")
    
    #check if exclude plastid==TRUE >call f_nuclid
    if (exclude_plastid==TRUE) {

      #call f_nuclid
      ls_geneid<-f_nuclid(querytext = querytext)
    }
    
    if (exclude_plastid==FALSE) {
      message("Message: exclude_plastid is optional, when plastid genes are to be included this field can be left empty."
            )
      message("No genes ID's are excluded.")
      #call f_allid
      ls_geneid<-f_allid(querytext = querytext)
      
    }
  }
      

  #if exclude nuclear is called check for FALSE >pass on, TRUE >call f_plastid, or ESLE >raise error
  if(
    missing(exclude_plastid) & 
    !missing(exclude_nuclear) & 
    missing(specify)){
    #test logic input
    f_logictest(logic_var = exclude_nuclear,name_var = "exclude nuclear")
    
    #check if exclude nuclear is TRUE
    if (exclude_nuclear==TRUE) {
      #call f_plastid function
      ls_geneid<-f_plastid(querytext = querytext)
    }
    if (exclude_nuclear==FALSE) {
      message("Message: exclude_nuclear is optional, when plastid genes are to be included this field can be left empty."
      )
      message("No genes ID's are excluded.")
      #call f_allid
      ls_geneid<-f_allid(querytext = querytext)
      
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
  
  
  }

  #if exclude nuclear | exclude plastid is TRUE and specify is specified >raise error
 
  
  if (
    (!missing(specify) &
    !missing(exclude_plastid)) |
    (!missing(specify) &
    !missing(exclude_nuclear))
    ){
    stop("Error: either use one of the exclude variables or use specify. Both is not supported at the moment.")
  }

#   #if both exclude nuclear and exlcude plastid is called
  if (
    
     !missing(exclude_plastid)&
     !missing(exclude_nuclear)
  ){ 
    #check for logic
    ##excl nucl
    f_logictest(logic_var = exclude_nuclear, name_var = "exclude_nuclear")
    
    ##excl plast
    f_logictest(logic_var = exclude_plastid, name_var = "exclude_plastid")

    # exclude only plastid
   if (
      exclude_plastid==TRUE &
      exclude_nuclear==FALSE
    ) {
      
      #call f_nucl
      ls_geneid<-f_nuclid(querytext = querytext)
      
    }
    
    # exclude only nuclear
    if (
      exclude_plastid==FALSE &
      exclude_nuclear==TRUE
    ) {
      
      #call f_nucl
      ls_geneid<-f_plastid(querytext = querytext)
      
    }
    
    # exclude neither >message they can be left empty >call f_allid
    if (
      exclude_plastid==FALSE &
      exclude_nuclear==FALSE
    ) {
      
      #call f_all +warning
      message("All gene ID's are being extracted. Be aware exclude_plastid and exclude nuclear are optional and can be left empty.")
      ls_geneid<-f_allid(querytext = querytext)
    }
    # exclude both >stop no genes to extract
    
     if (
      exclude_plastid==TRUE &
      exclude_nuclear==TRUE
    ) {
      
      #raise error >no genes to extract
      stop("Error: There are no ID's to extract when excluding both nuclear and plastid ID's.")
      
    }
    }
  

  
return(ls_geneid)

}
