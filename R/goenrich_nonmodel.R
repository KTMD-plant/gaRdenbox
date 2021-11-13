#goseq for non model species

#Author: Kilian Duijts



non_model_goseq<-function(  #fuction runs goseq for non-modelspecies
  genelist,                 #provide query list
  padj,                     #choose padj thresshold
  geneid2GO,                #provide named list. names=geneid, list elements goterms
  cdna                      #fasta file of sequences
){
  
  
  #get the length of each transcript
  bias_cdna<-nchar(cdna)
  
  
  #get the corresponding geneid
  names(bias_cdna)<-tibble(Sp=names(cdna)) %>%
    separate(Sp,into="Sp2",sep=" ",extra="drop") %>%
    dplyr::select(Sp2) %>% 
    as_vector()
  

  
  
  
  #get genelist to test
  tgl<-as.integer(
    names(bias_cdna) %in% genelist)
  
  names(tgl)<-names(bias_cdna)
  
  
  #calculate probability weighting
  pw<-nullp(tgl, bias.data = bias_cdna,plot.fit = TRUE)
  
  
  goseq_out<-goseq(pwf = pw, gene2cat = geneid2go)
  
  
  
  goseq_out[, "padj_BH"] <- p.adjust(goseq_out[, 2], method = "BH")
  
  GOenriched <- goseq_out[which(goseq_out[, "padj_BH"] < padj), ]
  y <- genelist[genelist %in% names(geneid2go)]
  
  if(nrow(GOenriched) != 0){
    GOenriched["Total"] <- nrow(pw)
    GOenriched["ClusterSize"] <- length(y)
    
    
    go2glgenes<-data.frame("category"=character(),
                           "geneid"=character())
    
    for (j in GOenriched$category) {
      
      a<-go2geneid[[j]]
      gogenes<-a[a%in%genelist]%>%
        unlist()
      
      
      df<-data.frame("category"=j,
                     "geneid"=paste(gogenes, 
                                    collapse = ","
                     )
      )
      
      go2glgenes<-rbind(go2glgenes, df)
      
    }
    
    GOenriched<-left_join(GOenriched, go2glgenes)
    
    
    
    
  }
  
  
  return(GOenriched)
  
}           

