# helpers.R
# helper functions for hapmap.R
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# SAIC-Frederick, Inc
# Created August 28, 2013
# Last Modified October 31, 2013

get.phased <- function(chr, pop, clean.duos = FALSE, dataDir = '../Data/')
{
    # read in and merge phased data (some individuals did not pass QC,
    #                                hence those that were duos or unrelated)

    # trios
    phased1 <- try(read.table(paste(dataDir, toupper(pop), '/phased/hapmap3_r2_b36_fwd.consensus.qc.poly.chr',
                                    chr, '_', tolower(pop), '.unr.phased', sep = ''),
                              header = TRUE, stringsAsFactors = FALSE, na.strings = '-'))
    # duos
    phased2 <- try(read.table(paste(dataDir, toupper(pop), '/phased/hapmap3_r2_b36_fwd.consensus.qc.poly.chr',
                                    chr, '_', tolower(pop), '.D.phased', sep = ''),
                              header = TRUE, stringsAsFactors = FALSE, na.strings = '-'))
    if(class(phased2) == 'data.frame' & clean.duos)
    {
        todrop <- c(grep('_NA', names(phased2), fixed = TRUE), grep('_0_A', names(phased2), fixed = TRUE))

        if(length(todrop) > 0)
            phased2 <- phased2[,-todrop]
    }
    # unrelated
    phased3 <- try(read.table(paste(dataDir, toupper(pop), '/phased/hapmap3_r2_b36_fwd.consensus.qc.poly.chr',
                                    chr, '_', tolower(pop), '.phased', sep = ''),
                              header = TRUE, stringsAsFactors = FALSE, na.strings = '-'))


    # merge...not all of these are defined in all populations :P
    if(class(phased1) == 'data.frame' & class(phased2) == 'data.frame' & class(phased3) == 'data.frame')
    {
        phased <- merge(merge(phased3, phased2), phased1)

        if(dim(phased)[1] != mean(dim(phased1)[1], dim(phased2)[1], dim(phased3)[1]))
            stop("Differing number of rows in chromosome ", chr)
    }else{
        if(class(phased1) == 'try-error')
        {
            if(class(phased2) == 'data.frame' & class(phased3) == 'data.frame')
            {
                phased <- merge(phased2, phased3)
                if(dim(phased)[1] != mean(dim(phased2)[1], dim(phased3)[1]))
                    stop("Differing number of rows in chromosome ", chr)
            }else{
                if(class(phased2) == 'try-error')
                    phased <- phased3
                if(class(phased3) == 'try-error')
                    phased <- phased2
            }
        }

        if(class(phased2) == 'try-error')
        {
            if(class(phased1) == 'data.frame' & class(phased3) == 'data.frame')
            {
                phased <- merge(phased1, phased3)
                if(dim(phased)[1] != mean(dim(phased1)[1], dim(phased3)[1]))
                    stop("Differing number of rows in chromosome ", chr)
            }else{
                if(class(phased1) == 'try-error')
                    phased <- phased3
                if(class(phased3) == 'try-error')
                    phased <- phased1
            }
        }

        if(class(phased3) == 'try-error')
        {
            if(class(phased2) == 'data.frame' & class(phased1) == 'data.frame')
            {
                phased <- merge(phased2, phased1)
                if(dim(phased)[1] != mean(dim(phased1)[1], dim(phased2)[1]))
                    stop("Differing number of rows in chromosome ", chr)
            }else{
                if(class(phased2) == 'try-error')
                    phased <- phased1
                if(class(phased1) == 'try-error')
                    phased <- phased2
            }
        }
    }

    # clean up "AFFY" junk at the beginning of the string
    phased$rsID[substr(phased$rsID, 1, 2) == 'AF'] <-
        unlist(sapply(strsplit(phased$rsID[substr(phased$rsID, 1, 2) == 'AF'], '__', fixed = TRUE), `[`, 2))

    attributes(phased)$chr <- chr
    attributes(phased)$pop <- pop

    return(phased)
}

get.LD <- function(chr, pop, keep)
{
    # read in LD data
    LD <- read.table(paste(dataDir, toupper(pop), '/LD/ld_chr', chr, '_', toupper(pop), '.txt', sep = ''),
                     stringsAsFactors = FALSE)
    names(LD) <- c('pos1', 'pos2', 'pop', 'rs1', 'rs2', 'Dprime', 'r2', 'lod', 'fbin')

    # only include those that are significant at 0.05 level and are in the phased data
    LD <- subset(LD, lod >= 3 & rs1 %in% keep & rs2 %in% keep)

    # clean up "AFFY" junk at the beginning of the string
    LD$rs1[substr(LD$rs1, 1, 2) == 'AF'] <-
        unlist(sapply(strsplit(LD$rs1[substr(LD$rs1, 1, 2) == 'AF'], '__', fixed = TRUE), `[`, 2))

    LD$rs2[substr(LD$rs2, 1, 2) == 'AF'] <-
        unlist(sapply(strsplit(LD$rs2[substr(LD$rs2, 1, 2) == 'AF'], '__', fixed = TRUE), `[`, 2))

    attributes(LD)$chr <- chr
    attributes(LD)$pop <- pop

    return(LD)
}
