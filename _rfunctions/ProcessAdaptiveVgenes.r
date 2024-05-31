#Function modified from ProcessAdaptiveFile.R written as part of GIANA. 
#Accessed from GIANA Github 17 September 2023 
# https://github.com/s175573/GIANA/blob/master/ProcessAdaptiveFile.R
ProcessAdaptiveVgenes <- function(input){
  
  #AAG 17 September 2023 
  Vgene_o=c(1,2,9,13:19,26:28,30)
  Vgene_o=as.character(Vgene_o)
  
  gsub('TCRBV[0]{0,1}','TRBV', input)->tmpV
  ## Multiple calls
  vv.m=grep('/',tmpV)
  if(length(vv.m)>0){
    TMP=strsplit(tmpV[vv.m],'/')
    tmpV[vv.m]=unlist(sapply(TMP,function(x)x[1]))
  }
  vv.m=grep(',',tmpV)
  if(length(vv.m)>0){
    TMP=strsplit(tmpV[vv.m],',')
    tmpV[vv.m]=unlist(sapply(TMP,function(x)x[1]))
  }
  tmpV=gsub('-0','-',tmpV)
  v_digit=grep('\\*',tmpV)
  if(length(v_digit)>0)tmpV[- v_digit] = paste(tmpV[- v_digit],'*01',sep='') else tmpV=paste(tmpV,'*01',sep='')
  ## 1. Orphan V genes do not have "-1", need to remove
  Vnumbers=gsub('TRBV','',tmpV)
  Vnumbers=gsub('\\*.+','',Vnumbers)
  Vnumbers=gsub('-.+','',Vnumbers)
  vv.o=which(Vnumbers %in% Vgene_o)
  tmpV1=tmpV
  tmpV1[vv.o]=gsub('-1','',tmpV1[vv.o])
  ## 2. Non-orphan V genes but without "-1", need to add
  vv.no=which(! Vnumbers %in% Vgene_o)
  vv.non=grep('-',tmpV1[vv.no])
  if(length(vv.non)>0)tmpV1[vv.no][-vv.non]=gsub('\\*01','-1*01',tmpV1[vv.no][-vv.non])
  output=tmpV1
  return(output)
}

