variantlist <- function(){
    savedatapath <- "../data/"
    ## path for tsv files 
    path = "../variantcalling/"
    
    
    files <-  list.files(path=path,pattern=".tsv$")
    allf <- files
    alllist <- c()
    for(i in 1:length(allf)){
        tmp <- paste(path,allf[i],sep="")
        oner <- read.delim(tmp)
        subj <- gsub(".tsv","",allf[i])
        oner <- cbind(oner,subj)
        cols <- colnames(oner)
        colsub <- c(which(grepl(paste(subj,".GT",sep=""),cols)),which(grepl(paste(subj,".AD",sep=""),cols)),max(which(grepl(subj,cols))), dim(oner)[2])
        colnames(oner)[colsub] <- c("GT","AD","Subject_INFO","Subject_ID")
        alllist <- rbind(alllist,oner)
    }
    save(alllist,file=paste(savedatapath,"alllist",sep=""))
    
    allV <- paste(alllist[,1],alllist[,2],alllist[,4],alllist[,5],sep="_")
    varT <- table(allV)
    save(varT,file=paste(savedatapath,"varT",sep=""))

}

filtering <- function(onelist){

    filters <- c("filtered","ExACfreq","VCFPASS","noneSegmentalDup","meta-SVM_PP2")
    filS <- matrix(FALSE,dim(onelist)[1],length(filters))
    colnames(filS) <- filters
    onelist <- cbind(onelist,filS)
    
    ## rare variants
    Ecut=0.005;
    onelist[is.na(onelist[,"AlleleFrequency.ExAC"]),"AlleleFrequency.ExAC"] <- 0
    onelist[onelist[,"AlleleFrequency.ExAC"]==".","AlleleFrequency.ExAC"] <- 0
    onelist[as.numeric(onelist[,"AlleleFrequency.ExAC"])< Ecut,"ExACfreq"] <- TRUE
    
    ### filtered more details
    badvars <- c('QD_Bad_SNP','FS_Bad_SNP','FS_Mid_SNP;QD_Mid_SNP','LowQuality','LowQD_Indel','LowQuality','VQSRTrancheSNP99.90to100.00','VQSRTrancheINDEL99.90to100.00')
    subs1 <- sapply(1:dim(onelist)[1],function(i){
        tmp <- unlist(strsplit(onelist[i,"FILTER"],";"))
        a1 <- intersect(tmp,badvars)
        a2 <- grepl("FS_Mid_SNP;QD_Mid_SNP",onelist[i,"FILTER"])
        (length(a1) > 0) | a2
    })
    subs1 <- !subs1
    
    subs2 <- sapply(1:dim(onelist)[1], function(i) {
        if(onelist[i,"SegmentalDuplication"] == "none"){ TRUE;
        }else{ tmp <- unlist(strsplit(onelist[i,"SegmentalDuplication"],","))[1]
               as.numeric(unlist(strsplit(tmp,":"))[2]) < 0.95
        }
    })
    onelist[subs1,"VCFPASS"] <- TRUE
    onelist[subs2,"noneSegmentalDup"] <- TRUE
    
    ### missense predicted by meta-SVM and PP2
    mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
    onelist[,"meta-SVM_PP2"] <- TRUE
    subs2 <- onelist[,"VariantClass"] %in% mis
    subs1 <- onelist[,"PP2prediction"]=="D" | onelist[,"MetaSVM"]=="D"
    #subs1 <- onelist[,"MetaSVM"]=="D"
    subs3 <- nchar(onelist[,"REF"]) != nchar(onelist[,"ALT"]) ### indels
    subs2 <- subs2 & !subs3
    onelist[subs2 & !subs1,"meta-SVM_PP2"] <- FALSE
    
    onelist[,"filtered"] <- onelist[,"ExACfreq"] & onelist[,"VCFPASS"] & onelist[,"noneSegmentalDup"] & onelist[,"meta-SVM_PP2"]
    
    onelist

}

transmmited <- function(onetrio,varlist,cangene){
    ## onetrio: at least contains one proband and one parents
    prolist <- varlist[varlist[,"Subject_ID"] %in% onetrio[2],]
    parlist <- varlist[varlist[,"Subject_ID"] %in% c(onetrio[3],onetrio[4]),]
    
    proV <- paste(prolist[,1],prolist[,2],prolist[,4],prolist[,5],sep="_")
    parV <- paste(parlist[,1],parlist[,2],parlist[,4],parlist[,5],sep="_")
    
    a <- length( unique( proV[(proV %in% parV) & (prolist[,"Gene"] %in% cangene)] ) )
    b <- length( unique( proV[(proV %in% parV) & !(prolist[,"Gene"] %in% cangene)] ) )
    c <- length( unique( proV[ !(proV %in% parV) & (prolist[,"Gene"] %in% cangene)] ) )
    d <- length( unique( proV[ !(proV %in% parV) & !(prolist[,"Gene"] %in% cangene)] ) )
    p <- fisher.test(matrix(c(a,b,c,d),2,2))$p.value
    
    c(a,b,c,d,p)
}
