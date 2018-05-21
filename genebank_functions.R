## ################### ##
## Gene Bank Curation  ##
## Narinder Singh      ##
## Updated 08/28/17    ##
## ################### ##

# function to compute error rate for regbs runs
ErrorRate <- function(reGBS.Run, compute.error.on) {
   cat('The error rate for', reGBS.Run, 'is',
       round(
          mean(
             colSums(1 - compute.error.on, na.rm = T) /
                colSums(!is.na(compute.error.on))) * 100, 2),
       '%')
}

alleleMatching <- function(allele.match=allele.match) {
    nS = ncol(allele.match)-9
    id = matrix(NA, nrow=nS, ncol=nS)
    for (i in 1:nrow(id)){
      id_ct <- rep(NA, length(i:nrow(id)))
      id_pc <- rep(NA, length(i:nrow(id)))
        for (j in i:ncol(id)){
          line1 = as.character(allele.match[,i+9])
          line2 = as.character(allele.match[,j+9])
          shared = line1!="N" & line2!="N" & line1!="H" & line2!="H"
          common = line1[shared] == line2[shared]
          id_pc[j-i+1] <- sum(common)/sum(shared)
          id_ct[j-i+1] <- sum(shared)
        }
      id[i,i:ncol(id)] <- id_ct
      id[i:ncol(id),i] <- id_pc
    }
      rownames(id)=colnames(allele.match)[-c(1:9)]
      colnames(id)=colnames(allele.match)[-c(1:9)]
      write.table(id, file="data/idFull.mat", quote=F, sep="\t")
      cat('Done computing. Identity matrix has been saved as a tab-delimited file named idFull.mat')
}

# function to create a data.frame with identical accessions
identicalAccessions <- function(id, rows = 1:nrow(id), cols = 1:ncol(id)) {

   for (i in 1:ncol(id)) { # loop to set replicates to 0
      id[which(rownames(id) %in% colnames(id)[i]), i] = 0
   }

   id[id > 1] = NA
   tmp.csv.out <- matrix(NA, nrow = length(rows), ncol = 2)
   tmpIdMat <- t(id[rows, cols])
   for (i in 1:nrow(tmpIdMat)) {
      tmp.csv.out[i, 1] <- rownames(tmpIdMat)[i]
      for (j in i:ncol(tmpIdMat)) {
         if (tmpIdMat[i, j] > 0.99)
            tmp.csv.out[i, 2] <- ifelse(is.na(tmp.csv.out[i, 2]),
                                        colnames(tmpIdMat)[j],
                                        paste(tmp.csv.out[i, 2], colnames(tmpIdMat)[j], sep=", "))
      }
   }
   return(tmp.csv.out)
}

# function to compute Nei's diversity index
gbs.diversity <- function(hap, het="H", missing="N") {

   a <- substr(hap$alleles,1,1)
   b <- substr(hap$alleles,3,3)

   # Recalculate allele counts
   missing <- rowSums(hap == missing, na.rm = T)
   het <- rowSums(hap == het, na.rm = T)
   alleleA <- rowSums(hap == a, na.rm = T)*2 + het
   alleleB <- rowSums(hap == b, na.rm = T)*2 + het

   # Minor allele frequency
   maf <- apply(cbind(alleleA, alleleB), 1, min) / apply(cbind(alleleA, alleleB), 1, sum)

   # Nei's diversity index
   cat("Computing Nei's diversity index...")
   nei <- mean(1-(maf^2)-(1-maf)^2, na.rm = T)
   cat(" Done")

   return(nei)

}
