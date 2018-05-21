## ################################ ##
## Gene Bank Curation	        	   ##
## Narinder Singh and Jesse Poland  ##
## Updated 05/19/17  	      	   ##
## #############################%## ##

# load functions
source('genebank_functions.R')

   # load data.table library for fast data reading
   library(data.table)

   #####################
   ## cluster generation
   #####################

   # read numericalized hap file
   hap01 <- fread("data/hap01.txt", header = T, data.table = F)

      # range of missing data
      range(hap01$missing)
      sum(hap01$missing < 0.5)

   hap01 <- hap01[hap01$missing < 0.5, ]
   hap01[1:5, 1:10]

   # keep only numericalized genotype columns
   geno <- hap01[,10:ncol(hap01)]
   geno <- t(geno)
   geno[1:5,1:5]

   # compute genetic distance
   ifelse(test = file.exists('data/distMat.RData'),
          yes = load('data/distMat.RData'), no = distMat <- dist(geno))

   # saving distance matrix because it's a computationally intensive step
   if (!file.exists('data/distMat.RData')) save(distMat, file = "data/distMat.RData")

   # cluster accessions
   hc <- hclust(distMat)

   # load ape library for tree construction
   library(ape)

   # convert hclust object to phylo object to be used with ape
   hc2 <- as.phylo(hc)

      # spliting accessions names by "_" delimiter to get simplified names
      idCols.simple <- sapply(strsplit(rownames(geno), "_"), function(x) x[1])
      sum(duplicated(idCols.simple))

      # compute total number of accessions in each collection
      wgrc <- unique(idCols.simple)[grep('TA', unique(idCols.simple))]; length(wgrc)
      pau <- unique(idCols.simple)[grep('PAU', unique(idCols.simple))]; length(pau)
      cimmyt <- unique(idCols.simple)[grep('GID', unique(idCols.simple))]; length(cimmyt)
      total.acc <- length(wgrc) + length(pau) + length(cimmyt); total.acc

   # genebank info
   genebank <- cbind(c(wgrc, pau, cimmyt), NA)
   genebank[, 2][grep('TA', genebank[, 1])] = 'WGRC'
   genebank[, 2][grep('PAU', genebank[, 1])] = 'PAU'
   genebank[, 2][grep('GID', genebank[, 1])] = 'CIMMYT'

   # cluster coloring
   edgecols <- cbind(NA, 1:nrow(hc2$edge), 'black')

      ## branch coloring
      # fill in the label names using 2nd col of hc2$edge
      edgecols[, 1] <- hc2$tip.label[hc2$edge[, 2]]
      edgecols[, 1] <- sapply(strsplit(edgecols[, 1], "_"), function(x) x[1])

      edgecols <- as.matrix(merge(edgecols, genebank, by = 1, all.x = T))
      edgecols <- edgecols[order(as.numeric(edgecols[, 2])), ]

      # coloring diff genebanks with different colors
      edgecols[, 3][edgecols[, 4] == "WGRC"] = "red"
      edgecols[, 3][edgecols[, 4] == "PAU"] = "blue"
      edgecols[, 3][edgecols[, 4] == "CIMMYT"] = "green"

      # tips coloring
      tipcols <- hc2$tip.label
      tipcols <- as.matrix(merge(tipcols, edgecols, by = 1, sort = F))

      # plotting tree
      library(ape)
      pdf("cluster_tauschii_genebanks.pdf", width = 20, height = 20)
         plot(hc2, type = 'u', lab4ut = "axial", label.offset = 1, cex = 0.25,
              edge.color = edgecols[, 3], edge.width = 2, tip.color = tipcols[, 3])
         title(cex.main = 2.5, line = -1.5,
               main = substitute(paste("Cluster analysis of WGRC, PAU and CIMMYT genebanks' ",
                                       italic("Aegilops tauschii"), " collections using GBS")))
         legend(0, 200, legend = c('WGRC', 'PAU', 'CIMMYT'),
                text.col = c("red", "blue", "green"), col = c("red", "blue", "green"),
                lty = 1, lwd = 5, cex = 1.5)
         text(400, 200, labels = 'Possible outliers', cex = 1.5)

         ### use if want to draw fan type tree
         # legend(180, 210, legend = c('WGRC', 'PAU', 'CIMMYT'),
         #        text.col = c("red", "blue", "green"), col = c("red", "blue", "green"),
         #        lty = 1, lwd = 5, cex = 1.5)
         # text(35, 7, labels = 'Point of divergence', cex = 1.5)
         # points(0, 0, col = 'gold', cex = 2, pch = 16)
         # points(-28, 72, col = 'red', cex = 2, pch = 16)
         # legend(-260, 210, legend = c('Point of divergence', 'Outlier divergence'),
         #        text.col = c("gold", "red"), col = c("gold", "red"),
         #        pch = 16, cex = 1.5)
      dev.off()

#######################################
## compute ID matrix by allele matching
#######################################
      hap <- fread("data/hapFile.txt", header = T, check.names = F, data.table = F)
      alleleMatching(allele.match = hap) # this is computationally intensive step

#################
## load id matrix
#################

   id <- as.matrix(read.table('data/idFull.mat', header = T, check.names = F))
   dim(id)
   all(colnames(id) == rownames(id)) # check if colnames and rownames are same

#######################################
## all accessions that are ever regbs'd
#######################################
   regbs.acc.all <- unique(idCols.simple[duplicated(idCols.simple)])
   sum(duplicated(regbs.acc.all))

   # get indices for regbs'd accessions and its replicates
   regbs.acc.indices <- which(idCols.simple %in% regbs.acc.all)

   # checking if all diagonal elements are equal to 1
   all(diag(id) == 1)

   # setting diagonal to zero
   if(all(diag(id) == 1)) diag(id) = NA

   # matrix for regbs'd accessions
   regbs.acc.mat <- id[regbs.acc.indices, regbs.acc.indices]
   all(rownames(regbs.acc.mat) == colnames(regbs.acc.mat))

   # setting values >1 to NA
   regbs.acc.mat[regbs.acc.mat > 1] = NA

   # getting simple names for regbs accessions
   names.regbs.mat <- sapply(strsplit(colnames(regbs.acc.mat), "_"), function(x) x[1])
   names.regbs.mat

   # setting regbs mat cols to simple names
   colnames(regbs.acc.mat) <- names.regbs.mat

   # retain only the original acc and remove replicate columns
   regbs.acc.mat <- regbs.acc.mat[, !duplicated(colnames(regbs.acc.mat))]

   # loop to set non-replicates to NA
      for (i in 1:ncol(regbs.acc.mat)) {
         regbs.indices <- which(!names.regbs.mat %in% colnames(regbs.acc.mat)[i])
         regbs.acc.mat[regbs.indices, i] = NA
      }

   # removing extra columns and rows for easy reading
   regbs.acc.mat <- regbs.acc.mat[rowSums(!is.na(regbs.acc.mat)) > 0,
                                  colSums(!is.na(regbs.acc.mat)) > 0]


###################
# regbs1 accessions
###################
   regbs1.acc <- colnames(id)[grep('regbs', colnames(id), ignore.case = T)]
   regbs1.acc <- regbs1.acc[-grep('regbs2', regbs1.acc, ignore.case = T)]
   regbs1.acc <- regbs1.acc[-grep('regbs3', regbs1.acc, ignore.case = T)]
   regbs1.acc.simple <- unique(sapply(strsplit(regbs1.acc, "_"), function(x) x[1]))

   # regbs1 matrix
   regbs1.mat <- regbs.acc.mat[rownames(regbs.acc.mat) %in% regbs1.acc,
                               colnames(regbs.acc.mat) %in% regbs1.acc.simple]
   rownames(regbs1.mat)
   colnames(regbs1.mat)

   # error rate for reGBS1
   ErrorRate('reGBS1', regbs1.mat)


#################################################################
## create a data.frame with accessions with lower id coefficients
#################################################################

   # set threshold for identity coefficient
   id.thresh = 0.99

   # variables to hold low id accession names and id coefficients
   low.match.acc <- NULL
   low.match.id <- NULL

   # loop to retrieve accession names and id coefficients
      for (i in 1:nrow(regbs1.mat)) {
         low.match.acc <- c(low.match.acc,
                            colnames(regbs1.mat)[which(regbs1.mat[i, ] < id.thresh)])
         low.match.id <- c(low.match.id,
                           regbs1.mat[i, ][which(regbs1.mat[i, ] < id.thresh)])
      }

   # combining variables to create a data frame
   low.match.acc.df <- data.frame(accession = as.character(low.match.acc),
                                  identity = as.numeric(low.match.id))
   # removing duplicated accessions
   low.match.acc.df <- low.match.acc.df[!duplicated(low.match.acc.df$accession), ]

   low.match.acc.df # print out the accessions and id coefficients


##########################################################################################
## create an object with individual matrices for low id accessions and regbs2 coefficients
##########################################################################################

   # initiate empty list
   low.id.mat.list <- list()

   # loop to retrieve matching reGBS accessions
      for (i in 1:nrow(low.match.acc.df)) {
         name <- as.character(low.match.acc.df$accession[i])
         mat <- id[grep(name, rownames(id)), grep(name, colnames(id))]
         # mat <- mat[grep('regbs3', rownames(mat), ignore.case = T, invert = T),
         #            grep('regbs3', colnames(mat), ignore.case = T, invert = T)]
         low.id.mat.list[[name]] <- mat
      }

   # print the list
   low.id.mat.list

   # list of low id samples to remove
   rem.low.acc <- c('TA10106', 'TA10108', 'TA1581_ReGBS', 'TA1589_ReGBS',
                    'TA2431_ReGBS', 'TA2468', 'TA2488_ReGBS')

#########################################
## original, regbs1 and regbs2 accessions
#########################################

   regbs2.acc <- colnames(id)[grep('regbs', colnames(id), ignore.case = T)]
   regbs2.acc <- regbs2.acc[-grep('regbs3', regbs2.acc, ignore.case = T)]
   regbs2.acc.simple <- unique(sapply(strsplit(regbs2.acc, "_"), function(x) x[1]))

   # regbs2 matrix
   regbs2.mat <- id[rownames(id) %in% regbs1.acc | rownames(id) %in% regbs2.acc,
                    colnames(id) %in% regbs1.acc | colnames(id) %in% regbs1.acc.simple]
   rownames(regbs2.mat)
   colnames(regbs2.mat)

   # setting non-replicates to NA
   name.cols <- sapply(strsplit(colnames(regbs2.mat), "_"), function(x) x[1])
   name.rows <- sapply(strsplit(rownames(regbs2.mat), "_"), function(x) x[1])
   for (i in 1:ncol(regbs2.mat)) {
      regbs.indices <- which(!name.rows %in% name.cols[i])
      regbs2.mat[regbs.indices, i] = NA
   }

   regbs2.mat <- regbs2.mat[!rownames(regbs2.mat) %in% rem.low.acc,
                            !colnames(regbs2.mat) %in% rem.low.acc]

   # setting values greater than zero to NA
   # setting identified very low coefficients (<0.85) to NA
   regbs2.mat[regbs2.mat > 1 | regbs2.mat < 0.85] = NA

# removing extra columns and rows for easy reading
   regbs2.mat <- regbs2.mat[rowSums(!is.na(regbs2.mat)) > 0,
                            colSums(!is.na(regbs2.mat)) > 0]

# error rate for regbs2
   ErrorRate('reGBS2', regbs2.mat)


# detecting heterogeneosity in TA1714
   TA1714 <- grep('1714_regbs3', colnames(id), ignore.case = T)
   TA2457 <- grep('2457_regbs', colnames(id), ignore.case = T)

   id[TA1714, TA1714]
   id[TA2457, TA2457]

# overall matching accessions for all three genebanks
   id2 <- id

   id2.row.names <- sapply(strsplit(rownames(id2), "_"), function(x) x[1])
   id2.col.names <- sapply(strsplit(colnames(id2), "_"), function(x) x[1])
   all(id2.row.names == id2.col.names)

   rownames(id2) <- id2.row.names
   colnames(id2) <- id2.col.names

   csv.out <- identicalAccessions(id2)

# filtering the csv.out
   # create an empty list
   lst <- list()

   # loop to parse out dup accessions
   for (i in 1:nrow(csv.out)) {
     lst[i] <- strsplit(csv.out[i, 2], ', ')
   }

   # keeping only unique accessions per group
   lst <- lapply(lst, unique)

   # combining unique accessions
   lst.to.mat <- sapply(lst, function(x) paste(x, collapse = ', '))
   csv.out[, 2] <- lst.to.mat

   # getting list of accessions in the second colums, removing dups
   lst = unique(na.omit(unlist(lst)))

   # removing the rows if any accessions exists in first column
   csv.out <- csv.out[!csv.out[, 1] %in% lst, ]
   csv.out <- csv.out[!duplicated(csv.out[, 1]), ]

   # checking total number of accessions
   nrow(csv.out) + length(lst)

   # group accessions in one column
   grp.accessions <- apply(csv.out, 1, function(x) paste(x, collapse = ', '))
   grp.accessions <- gsub(', NA', '', grp.accessions)

   # sort the accessions groups
   grp.accessions <- grp.accessions[order(grp.accessions)]

   # compute group size
   grp.size <- sapply(strsplit(grp.accessions, ', '), length)
   grp.size.dist <- table(grp.size)

   # plot group size distribution
   pdf('barplot_grp_size.pdf', width = 11, height = 8.5)
   bplot <- barplot(grp.size.dist, ylim = c(0, max(grp.size.dist) + 50),
                    cex.axis = 1.3)
   text(bplot, grp.size.dist + 15, labels = grp.size.dist, col = 'red')
   mtext(text = 'Number of individuals in each group', side = 1, line = 2.5, cex = 1.5)
   mtext(text = 'Frequency', side = 2, line = 2.5, cex = 1.5)
   dev.off()

##########################################
## count group categories for venn diagram
##########################################
   grp.cnt <- matrix(NA, nrow = length(grp.accessions), ncol = 4)
      for (i in 1:nrow(grp.cnt)) {
         grp.cat = NULL
         split.grp <- strsplit(grp.accessions[i], ', ')

            wgrc.len <- length(grep('TA', split.grp[[1]]))
            grp.cnt[i, 1] = wgrc.len
               if (wgrc.len > 0) grp.cat <- 'W'

            pau.len <- length(grep('PAU', split.grp[[1]]))
            grp.cnt[i, 2] = pau.len
               if (pau.len > 0) grp.cat <- paste(grp.cat, 'P', sep = '')

            cimmyt.len <- length(grep('GID', split.grp[[1]]))
            grp.cnt[i, 3] = cimmyt.len
               if (cimmyt.len > 0) grp.cat <- paste(grp.cat, 'C', sep = '')

            grp.cnt[i, 4] = grp.cat
      }

   gc()

   # assign group numbers
   grp.num <- paste('Grp', 1:length(grp.accessions), sep = '')

   # combine the groups and size in a data frame
   grp.accessions.data <- as.data.frame(cbind(grp.num, grp.size, grp.cnt, grp.accessions))
   colnames(grp.accessions.data) <- c('Grp#', 'Grp_size', 'WGRC',
                                      'PAU', 'CIMMYT', 'Grp_cat', 'Matching_accessions')

   # check group size equal to total of each collection
   total <- as.numeric(levels(grp.accessions.data$Grp_size))[grp.accessions.data$Grp_size]
   ind.coll.count <- as.matrix(grp.accessions.data[, 3:5])
   computed.total <- apply(ind.coll.count, 1, function(x) sum(as.numeric(x)))

   all(total == computed.total)

   # write out the accession groupings
   write.csv(grp.accessions.data, file = 'Matching_accessions.csv', row.names = F, quote = T)

#################
## ploting graphs
#################

   # read in the accessions data
   grp.accessions.data <- read.csv(file = 'Matching_accessions.csv',
                                   header = T, as.is = T, check.names = F)

   # total of each collection accessions should be equal to total number of accessions
   total.acc == length(unique(idCols.simple))

   # compute number of groups within collections
   wgrc.grps <- sum(grp.accessions.data$WGRC > 0); wgrc.grps
   pau.grps <- sum(grp.accessions.data$PAU > 0); pau.grps
   cimmyt.grps <- sum(grp.accessions.data$CIMMYT > 0); cimmyt.grps

   # calculate percent uniqueness
   whole.uniq = (nrow(grp.accessions.data) / total.acc) * 100; whole.uniq
   wgrc.uniq = (wgrc.grps / length(wgrc)) * 100; wgrc.uniq
   pau.uniq = (pau.grps / length(pau)) * 100; pau.uniq
   cimmyt.uniq = (cimmyt.grps / length(cimmyt)) * 100; cimmyt.uniq

   #############################################################
   ## average pairwise identity (duplication) within collections
   #############################################################
   mean(c(100 - wgrc.uniq, 100 - pau.uniq, 100 - cimmyt.uniq))

   # barplot
   library(ggplot2)
   pdf('barplot_collections.pdf')

   alpha.val = 0.4
   uniqueness <- c(whole.uniq, wgrc.uniq, pau.uniq, cimmyt.uniq)
   bplot <- barplot(uniqueness, names.arg = c("Whole", "WGRC", "PAU", "CIMMYT"),
           col = alpha(c('black', 'red', 'blue', 'green'), alpha.val),
           ylim = c(0, 100), width = 10, xlim = c(0, 50), cex.names = 1.5, cex.axis = 1.5)
   text(bplot, uniqueness + 6, labels = round(uniqueness, 2), cex = 1.5)
   mtext("% Unique accessions", 2, cex = 1.5, line = 2.5)

   dev.off()

   # group frequency
   grp.cnt.tab <- as.data.frame(t(as.matrix(table(as.character(grp.accessions.data$Grp_cat)))))

   # venn diagram
   library(venneuler)
   vd <- venneuler(c(W = grp.cnt.tab$W, P = grp.cnt.tab$P, C = grp.cnt.tab$C,
                     "W&P" = grp.cnt.tab$WP, "W&C" = grp.cnt.tab$WC, "P&C" = grp.cnt.tab$PC,
                     "W&P&C" = grp.cnt.tab$WPC))
   vd$labels <- ''
   pdf('venn_collections.pdf')
      cex = 1.3
      plot(vd, col = c('red', 'blue', 'green'), alpha = alpha.val)

         text(x = 0.85, y = 0.15, cex = cex, labels = 'WGRC')
            text(x = 0.75, y = 0.25, cex = cex, labels = grp.cnt.tab$W)

         text(x = 0.18, y = 0.5, cex = cex, labels = 'PAU')
            text(x = 0.27, y = 0.5, cex = cex, labels = grp.cnt.tab$P)

         text(x = 0.8, y = 0.7, cex = cex, labels = 'CIMMYT')
            text(x = 0.55, y = 0.75, cex = cex, labels = grp.cnt.tab$C)

         text(x = 0.43, y = 0.5, cex = cex, labels = grp.cnt.tab$WPC)
         text(x = 0.38, y = 0.43, cex = cex, labels = grp.cnt.tab$WP)
         text(x = 0.37, y = 0.58, cex = cex, labels = grp.cnt.tab$PC)
         text(x = 0.6, y = 0.55, cex = cex, labels = grp.cnt.tab$WC)

   dev.off()

   #############
   ## histograms
   #############

   ### overall pairwise identity
   pdf('pIBS_histograms.pdf', width = 11, height = 8.5)
   par(mfrow=c(2,2))

   hist(id[lower.tri(id)], xlab = 'pIBS',
        main = 'pIBS distribution for three genebanks collectively')
   abline(v = mean(id[lower.tri(id)]), col = 'red', lty = 2)

   wgrc.hist <- id[rownames(id) %in% wgrc, colnames(id) %in% wgrc]
   hist(wgrc.hist[lower.tri(wgrc.hist)], xlab = 'pIBS',
        main = 'pIBS distribution for WGRC genebank')
   abline(v = mean(wgrc.hist[lower.tri(wgrc.hist)]), col = 'red', lty = 2)

   cimmyt.hist <- id[rownames(id) %in% cimmyt, colnames(id) %in% cimmyt]
   hist(cimmyt.hist[lower.tri(cimmyt.hist)], xlab = 'pIBS',
        main = 'pIBS distribution for CIMMYT genebank')
   abline(v = mean(cimmyt.hist[lower.tri(cimmyt.hist)]), col = 'red', lty = 2)

   pau.hist <- id[rownames(id) %in% pau, colnames(id) %in% pau]
   hist(pau.hist[lower.tri(pau.hist)], xlab = 'pIBS',
        main = 'pIBS distribution for PAU genebank')
   abline(v = mean(pau.hist[lower.tri(pau.hist)]), col = 'red', lty = 2)

   dev.off()

##################
## diversity stats
##################
   # load hapfile generated from general code
   hap = fread("data/hapFile.txt", header = T, data.table = F)
   hap[1:5,1:10]
   colnames(hap)

   wgrc.names <- colnames(hap)[grep('TA', colnames(hap))]
   pau.names <- colnames(hap)[grep('PAU', colnames(hap))]
   cimmyt.names <- colnames(hap)[grep('GID', colnames(hap))]

   hap.wgrc <- hap[, c(1:9, which(colnames(hap) %in% wgrc.names))]
   hap.pau <- hap[, c(1:9, which(colnames(hap) %in% pau.names))]
   hap.cimmyt <- hap[, c(1:9, which(colnames(hap) %in% cimmyt.names))]

   nei.wgrc <- gbs.diversity(hap = hap.wgrc)
   nei.pau <- gbs.diversity(hap = hap.pau)
   nei.cimmyt <- gbs.diversity(hap = hap.cimmyt)

   # print out nei
   nei.wgrc; nei.pau; nei.cimmyt

######
## END
######
