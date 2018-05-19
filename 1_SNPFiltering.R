## ################## ##
## SNP filtering      ##
## Narinder Singh	  ##
## Revised 4/28/2017  ##
## ################## ##

# read all 3 hap files
hap1 <- read.table(file = "data/_170419_tauschii_wgrc1_20170429.hap", header = T, check.names = F, stringsAsFactors = F)
hap2 <- read.table(file = "data/_170419_tauschii_wgrc2_20170429.hap", header = T, check.names = F, stringsAsFactors = F)
hap3 <- read.table(file = "data/_170419_tauschii_wgrc3_20170429.hap", header = T, check.names = F, stringsAsFactors = F)

# check if each hapfile has same number of columns (should be TRUE)
(ncol(hap1) == ncol(hap2)) == (ncol(hap2) == ncol(hap3))

# number of columns
ncol(hap1)

# replacing $center column of hap1,hap2,hap3 with 1,2,3, respectively to keep them in order
hap1$center = 1; hap2$center = 2; hap3$center = 3

# checking if columns of all hap files are same (should be TRUE)
all(colnames(hap1) == colnames(hap2)) == all(colnames(hap1) == colnames(hap3))

# create new object with all hap files combined by rows
joined <- rbind(hap3, hap2, hap1)

# check number of rows in the joined and total number of rows in hapfiles (should be TRUE)
nrow(joined) == sum(nrow(hap1), nrow(hap2), nrow(hap3))

# free some memory
rm(hap1, hap2, hap3)
gc()

# create an onject with only rs and assembly columns
rs_pos <- joined[,c(1,6)]
rs_pos[1:5,]

# finding duplicated tag and assembly combination
dup <- duplicated(rs_pos)
sum(dup)

# removing duplicated entries
noDup <- joined[!dup,]
noDup[1:5,1:15]

# check number of non-duplicated tags
nrow(noDup)

# arranging rs and assembly in order
odr <- order(as.vector(noDup$rs), as.vector(noDup$assembly))

# arranging rows in according to the order of noDup
hap <- noDup[odr, ]
hap[1:5,1:15]
dim(hap)

# replacing column names
colnames(hap)[colnames(hap) == "assembly"] = "snp_pos"
colnames(hap)[colnames(hap) == "protLSID"] = "alleleA"
colnames(hap)[colnames(hap) == "assayLSID"] = "alleleB"
colnames(hap)[colnames(hap) == "panelLSID"] = "het"
colnames(hap)[colnames(hap) == "QCcode"] = "missing"
hap[1:5,1:20]

# check the blank wells
blank <- hap[, grepl("blank", colnames(hap), ignore.case = T)]
n.cols <- ncol(blank)

if (is.null(n.cols)) {
   blank <- as.character(blank)
   if(is.null(n.cols)) cat('Blank well has', sum(blank != "N"), 'sequence reads')
   } else if (n.cols == 0) {
         'There are no blank wells'
      } else {
         blankSum <- apply(blank!="N", 2, sum)
         if (exists('blankSum')) print(as.matrix(blankSum))
}

# removing specific columns
hap <- hap[,colnames(hap) != "chrom"]
hap <- hap[,colnames(hap) != "pos"]
hap <- hap[,colnames(hap) != "center"]
hap <- hap[,colnames(hap) != "strand"]
hap <- hap[, !grepl("blank", colnames(hap), ignore.case = T)]
hap[1:5, 1:15]
dim(hap)

# getting alleles
a <- substring(hap$alleles, 1, 1)
b <- substring(hap$alleles, 3, 3)
sum(a == b) # should be zero

# computing correct allele counts
alleleA <- rowSums(hap == a)
alleleB <- rowSums(hap == b)
het <- rowSums(hap == 'H')
missing <- rowSums(hap == 'N') / (ncol(hap) - 8)

# replacing hap columns with correct counts
hap$alleleA <- alleleA
hap$alleleB <- alleleB
hap$het <- het
hap$missing <- missing
hap[1:5, 1:15]

# calculating minor allele freq and percent HET
MAF <- apply(cbind((2*apply(cbind(hap$alleleA, hap$alleleB), 1, min)), hap$het), 1, sum) /
            (2*apply(cbind(hap$alleleA, hap$alleleB, hap$het), 1, sum))
HET <- hap$het / apply(cbind(hap$alleleA, hap$alleleB, hap$het), 1, sum)

hap <- cbind(hap[,c(1:7)], MAF, HET, hap[, c(8:ncol(hap))])
hap[1:5, 1:15]

# writing hap file in text file format
write.table(hap, file = "data/hapFile.txt" , col.names = T, row.names = F, sep = "\t", quote = F)

# remove everything
rm(list = ls())

# read in the final hapfile
hap <- read.table(file = "data/hapFile.txt", header = T, check.names = F, stringsAsFactors = F)
hap[1:5, 1:10]
dim(hap)

# basic histograms
pdf(file = "basic_histograms.pdf")

# histogram of minor allele frequency
hist(hap$MAF, main = "Minor allele frequency", xlab = "MAF", ylab = "#SNPs")

# histogram of HET SNPs
hist(hap$HET, main = "HET SNPs frequency", xlab = "% heterozygosity", ylab = "#SNPs")

# histogram of missing SNP data
hist(as.numeric(hap$missing), main = "Missing SNP data", xlab = ".prop absent", ylab = "#SNPs")

# histogram of Lines data for #SNPs
lineData <- colSums(hap[, -c(1:9)] != "N") / nrow(hap)
hist(lineData, main = "Line data", xlab = ".prop SNPs present", ylab = "# Lines")

dev.off()

# removing extra arguments to save space
rm(lineData)

# convert to alleles to numeric calls
hap01 <- hap
hap01[, 10:ncol(hap01)] = NA
hap01[1:5, 1:50]

# getting allele states
a = substring(hap$alleles, 1, 1)
b = substring(hap$alleles, 3, 3)

# changing first allele to second if second is major
a[hap$alleleA < hap$alleleB] = substring(hap$alleles, 3, 3)[hap$alleleA < hap$alleleB]
b[hap$alleleA < hap$alleleB] = substring(hap$alleles, 1, 1)[hap$alleleA < hap$alleleB]
sum(a == b)

# convert major allele to -1 and minor allele to 1, hets to 0
hap01[hap == a] = -1
hap01[hap == b] = 1
hap01[hap == "H"] = 0
hap01[1:10, 1:20]

# remove a, b
rm(a, b)

# writing hap01 file into a text format
write.table(hap01, file = "data/hap01.txt" , col.names = T, row.names = F, sep = "\t", quote = F)

# read hap01.txt to see if looks correct
hap01 = read.table(file = "data/hap01.txt", header = T, check.names = F, stringsAsFactors = F)
hap01[1:5, 1:10]
dim(hap01)

############################
##  End of SNP filtering  ##
############################
