#goseq for non model species

#Author: Kilian Duijts


#' @title Gene Ontology testing for non-model species
#' @description 
#' Re-application of the goseq ("https://bioconductor.org/packages/release/bioc/html/goseq.html") package for use on non model species.
#' This is based on a tutorial ("https://knozue.github.io/post/2018/10/17/over-representation-analysis-3-goseq-with-non-model-go-term.html") by knozue
#' @details 
#' The function makes it easier to test gene sets of non-model species.
#' 
#' It takes the following parameters:
#' 
#' @param querylist A list with for each element a non-model species gene id.
#' @param padj Numerical value specifying the cut-off for the adjusted pvalue
#' 
#' # method <currently not available>
#' 
#' @param geneid2GO A named list for which the names are geneids and the list elements consist of GO-terms (like so: GO:0006355). Please make sure your geneid's are identical to the query genes.
#' @param cDNA fasta file containing the non-model species geneid's and their sequences. The gene lenght is used for the probabily weight; normalization for gene length.
#' 
#' 



non_model_goseq<-function(  #fuction runs goseq for non-modelspecies
  querylist,                 #provide query list
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
  

  
  
  
  #get querylist to test
  tgl<-as.integer(
    names(bias_cdna) %in% querylist)
  
  names(tgl)<-names(bias_cdna)
  
  
  #calculate probability weighting
  pw<-goseq::nullp(tgl, bias.data = bias_cdna,plot.fit = TRUE)
  
  
  goseq_out<-goseq::goseq(pwf = pw, gene2cat = geneid2go)
  
  
  
  goseq_out[, "padj_BH"] <- p.adjust(goseq_out[, 2], method = "BH")
  
  GOenriched <- goseq_out[which(goseq_out[, "padj_BH"] < padj), ]
  y <- querylist[querylist %in% names(geneid2go)]
  
  if(nrow(GOenriched) != 0){
    GOenriched["Total"] <- nrow(pw)
    GOenriched["ClusterSize"] <- length(y)
    
    
    go2glgenes<-data.frame("category"=character(),
                           "geneid"=character())
    
    for (j in GOenriched$category) {
      
      a<-go2geneid[[j]]
      gogenes<-a[a%in%querylist]%>%
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

