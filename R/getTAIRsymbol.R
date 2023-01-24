# Functions to swap between agi codes



#function to go from agi code to symbol

##gets gene sybol from atgeneid input
f_get_at_symb<-function(atgeneid){
  
  require("org.At.tair.db")
  
  #at.org can not deal with the .1 or .2 splicing variants in the genid;
  ##Here I use the gettairid function from the gaRdenplot package to extract just 1 geneid
  ids<-data.frame("geneid"=gaRdenbox::gettairid(atgeneid))
  
 
    
    #query AT.org
    ids$map<-AnnotationDbi::mapIds(org.At.tair.db::org.At.tair.db,
                                   keys = ids$geneid,
                                   column = c('SYMBOL'), keytype = 'TAIR')
    
    #create one consistent column; replace NA symbols with agi code
    ids$map[is.na(ids$map)]<-ids$geneid[is.na(ids$map)]
    
    
    #consistent type case
    ids$map<-toupper(ids$map)
    
    message("Arabidopsis geneid found, converting to symbol")
    
  
  
  return(ids$map)
  
}

#function to get geneid from symbols; still needs to be implemented for the input.

f_get_atgenid_from_symb<-function(atsymb){
  ids<-data.frame("symbol"=atsymb)
  
  
  #query AT.org
  ids$map<-AnnotationDbi::mapIds(org.At.tair.db::org.At.tair.db,
                                 keys = ids$atsymb,
                                 column = c('TAIR'), keytype = 'SYMBOL')
  
  #create one consistent column
  ids$map[is.na(ids$map)]<-ids$atsymb[is.na(ids$map)]
  
  
  #consistent type case
  ids$map<-toupper(ids$map)
  
  return(ids$map)
}


