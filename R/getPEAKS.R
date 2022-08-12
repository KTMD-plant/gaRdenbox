
#autor: Kilian Duijts
#date 20220812


#functions used in GETPEAKS

#helper function


f_process_peaks<-function(fl_name, colid, chromosome, start_position_gene, length_upstream){
  
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



#' @title get_peaks
#' 
#' @description 
#' Loops through .NarrowPeaks files in input folder and extracts the peaks in the specified (promotor) secition of the genome.
#' 
#' @details 
#' This function has been independently developed to query files from the O'Malley paper ( doi: 10.1016/j.cell.2016.04.038).
#' Download input data from http://neomorph.salk.edu/dap_web/pages/browse_table_aj.php; in .NarrowPeaks format
#' 
#' @param fl_folder directory containining the .NarrowPeak files. The function runs recursive by default (non-recursive currently not supported). When using this function on O'Malley's dataset set the path to "dap_download_may2016_peaks".
#' @param chromosome chromosome number of the section of interest.
#' @param start_position_gene The position on the chromosome used as starting position of the section of interest. If only interested in the promotor use startcodon. If also interested in potenial binding sites in the introns/coding regions provide the position of the stopcodong or 3' UTR.
#' @param length_upstream The length of promotor sequeance of interest. Provide a negative number when querying for a gene in reverse orientation.
#'   
#' @return dataframe containing peaks and associated transcriptionfactors
#'   
#' @example 
#'  
#'  #to be added
#'  
#' @author Kilian Duijts
#'  
#' @export
#' 
#' 
#' 



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
    df<-f_process_peaks(fl_name = i,
                  colid = colid,
                  chromosome = chromosome,
                  start_position_gene = start_position_gene,
                  length_upstream = length_upstream)
    
    df_total_peaks<-rbind(df_total_peaks, df)
  }
  
  
  df_total_peaks<-df_total_peaks%>%na.omit()
  
  
  #extract gene ids/synbols
  df_total_peaks$geneid<-stringr::str_extract(string = df_total_peaks$tf, pattern = "^.+(?=_col)")
  
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
  df_total_peaks<-dplyr::left_join(df_total_peaks, ids)
  return(df_total_peaks)
}


#' @title get_tf_from_peaks
#' 
#' @description extracts only the unique geneids and symbols of transcription factors binding in the specified section.
#' 
#' @details This function uses get_peaks but returns a dataframe with only the unique geneids and names for the peaks in the section of interest.
#' 
#' parameters are identical to get_peaks.
#' 
#' @param fl_folder directory containining the .NarrowPeak files. The function runs recursive by default (non-recursive currently not supported). When using this function on O'Malley's dataset set the path to "dap_download_may2016_peaks".
#' @param chromosome chromosome number of the section of interest.
#' @param start_position_gene The position on the chromosome used as starting position of the section of interest. If only interested in the promotor use startcodon. If also interested in potenial binding sites in the introns/coding regions provide the position of the stopcodong or 3' UTR.
#' @param length_upstream The length of promotor sequeance of interest. Provide a negative number when querying for a gene in reverse orientation.
#'   
#'   @return dataframe containing associated transcription factor geneid's and symbols.
#'   
#'   @author Kilian Duijts
#' 
#' 
#' @export

#get function including

get_tf_from_peaks<-function(fl_folder, chromosome, start_position_gene, length_upstream){
  
  pk<-get_peaks(fl_folder = fl_folder,
                chromosome = chromosome,
                start_position_gene = start_position_gene, 
                length_upstream = length_upstream)
  
  
  tf<-pk[, c( "geneid", "map")]%>%unique()
  
  
  
  return(tf)
}

