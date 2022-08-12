
#autor: Kilian Duijts
#date 20220812


#functions used in GETPEAKS




f_process<-function(fl_name, colid, chromosome, start_position_gene, length_upstream){
  
  #get tf name
  tfname<-dirname(fl_name)%>%dirname()%>%basename()
  
  
  
  
  #calculate promotor boundries
  
  promstart= start_position_gene-length_upstream
  #promstart= qpos-2000
  
  #load file
  df<-read.delim(file = fl_name,sep = "\t",header = F)
  
  #assign column names
  colnames(df)<-colid
  
  
  df$tf<-tfname
  df$promotordistance<-start_position_gene-df$chromend
  
  df_peaks<-df%>%filter(chrom==paste("chr", chromosome, sep=""))%>%
    filter(chromstart>=promstart &
             chromend<=start_position_gene
    )
  
  
  
  return(df_peaks)
}







get_peaks<-function(fl_folder, chromosome, start_position_gene, length_upstream){
  #list files
  
  fl<-list.files(path = fin,pattern = ".narrowPeak",
                 full.names = T,recursive = T)
  
  #header ids
  colid<-c("chrom", "chromstart", "chromend", "name", "score", "strand", "signalvalue", "pvalue", "qvalue","peak")
  
  
  
  #initalize df
  
  df_total_peaks<-as.data.frame(matrix(nrow = 1, ncol = length(colid)+2))
  
  colnames(df_total_peaks)<-c(colid, "tf", "promotordistance")
  
  
  #itterate through files
  for (i in fl) {
    df<-f_process(fl_name = i,
                  colid = colid,
                  chromosome = chromosome,
                  start_position_gene = start_position_gene,
                  length_upstream = length_upstream)
    
    df_total_peaks<-rbind(df_total_peaks, df)
  }
  
  
  df_total_peaks<-df_total_peaks%>%na.omit()
  
  
  #extract gene ids/synbols
  df_total_peaks$geneid<-str_extract(string = df_total_peaks$tf, pattern = "^.+(?=_col)")
  
  require(org.At.tair.db)
  
  
  #get unique subset
  ids<-data.frame("geneid"=df_total_peaks$geneid%>%unique())
  
  
  #query AT.org
  ids$map<-AnnotationDbi::mapIds(org.At.tair.db::org.At.tair.db,
                                 keys = ids$geneid,
                                 column = c('TAIR'), keytype = 'SYMBOL') 
  
  #create one consistent column 
  ids$map[is.na(ids$map)]<-ids$geneid[is.na(ids$map)]
  
  
  #consistent type case
  ids$map<-toupper(ids$map)
  
  #join back to pk
  df_total_peaks<-left_join(df_total_peaks, ids)
  return(df_total_peaks)
}

#get function including

get_tf_from_peaks<-function(fl_folder, chromosome, start_position_gene, length_upstream){
  
  pk<-get_peaks(fl_folder = fl_folder,
                chromosome = chromosome,
                start_position_gene = start_position_gene, 
                length_upstream = length_upstream)
  
  
  tf<-pk[, c( "geneid", "map")]%>%unique()
  
  
  
  return(tf)
}

