source("~/.Rprofile")
#### transmitted test
#### candiate gene list
lf1 <- read.delim("../data/Combined_CLCP_gene_list.txt",blank.lines.skip = TRUE)
lf2 <- read.delim("../data/Summary_CLCP.txt",blank.lines.skip = TRUE)
cangene <- union(lf1[,"Gene"],lf2[,"Gene"])
cangene <- setdiff(cangene,"")
qwt(cangene,file="../data/candidateG.txt")

#### pedigree information
pedf <- read.delim("../data/CLCP_VCF_Key.txt")
pedf <- pedf[order(pedf[,3]),]
peds <- matrix(-9,dim(pedf)[1],6)
peds[,1] <- pedf[,3]
peds[,2] <- pedf[,1]
peds[,5] <- "other"
aff_sta <- c("CLCP","CP","CL","CLCP (syndromic)","Submucosal CP","BCL","BCLCP")
un_sta <- c("UCL","UCL-L","UCL-R","UCLCP","UCP-R")
for(i in 1:dim(pedf)[1]){
    if(pedf[i,5]=="YES"){
        tmp <- pedf[pedf[,3]==peds[i,1],]
        if(any(tmp[,6]=="Father")){
            peds[i,3] <- tmp[tmp[,6]=="Father",1]
        }
        if(any(tmp[,6]=="Mother")){
            peds[i,4] <- tmp[tmp[,6]=="Mother",1]
        }
    }
}
peds[pedf[,7] %in% aff_sta ,6] <- 2
peds[pedf[,7] %in% un_sta ,6] <- 1
qwt(peds,file="../data/pedigreeALL.ped")
subs <- peds[,3]!=-9| peds[,4]!=-9
ids <- setdiff(union(peds[subs,2],union(peds[subs,3],peds[subs,4])),-9)
qwt(peds[peds[,2] %in% ids,],file="../data/pedigreeSUB.ped")

#### sample IDs
c <- readLines(con=file("../data/CLCP_June2015.vcf","r"),n=202)
samples <- unlist(strsplit(c[201],"\t"))
samples <- samples[10:61]
qwt(samples,file="../data/SamplesID.txt")


#### variant lists

source("misc.R")
load("../data/alllist")
varlist <- filtering(alllist)
varlist <- varlist[varlist[,"filtered"],]

onetrio <- peds[49,]
transmmited(onetrio,varlist,cangene)

