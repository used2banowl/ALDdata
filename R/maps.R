# maps.R
# Read in rutgers maps from Matise et al and HapMap recombination maps
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created November 24, 2006
# Last Modified October 29, 2013

rutgers36 <- NULL
rutgers37 <- NULL
hapmap36 <- NULL

compare.maps <- FALSE

for(i in 1:23)
{
    ### Rutgers build 36 ###
    tmp <- read.table(paste('../data/Maps/rutgers_map_b36/chr', i,'.map2', sep=''),
                      sep='\t', na.strings=c('', ' '), header=TRUE)

    tmp$chr <- i

    if(compare.maps)
    {
        plot(tmp[[6]], tmp[[7]], type = 'l', ylab = 'cM', xlab = 'bp')
    }

    if(i == 23)
        tmp$Sex.averaged_map_position <- tmp$Female_map_position

    if(is.null(rutgers36)){
        rutgers36 <- tmp
    }else{
        rutgers36 <- merge(rutgers36, tmp, all=TRUE)
    }

    ### Rutgers build 37 ###
    # this turns out to be exactly the same either with interpolated or with only the backbone
    # save on the space and interpolate this myself!
    system(paste('more ../data/Maps/rutgers_map_b37/RUMapv3_B137_chr', i, '.txt | grep ackbone',
                 ' > tmp.txt', sep = ''))
    tmp <- read.table('tmp.txt', sep = '\t', na.strings = c('', ' ', 'NA'), header = TRUE)

    tmp$chr <- i

    if(i == 23)
        tmp$Sex_averaged_start_map_position <- tmp$Female_start_map_position

    if(compare.maps)
    {
        lines(tmp[[6]], tmp[[7]], col = 'blue')
    }

    if(is.null(rutgers37)){
        rutgers37 <- tmp
    }else{
        rutgers37 <- merge(rutgers37, tmp, all=TRUE)
    }

    ### HapMap build 36 ###
    if(i < 23)
    {
        tmp <- read.table(paste('../data/Maps/HapMap/genetic_map_chr', i, '_b36.txt', sep = ''),
                          header = TRUE)

        if(compare.maps)
        {
            lines(tmp[[1]], tmp[[3]], col = 'green')
        }

        tmp[[2]] <- NULL
        tmp$chr <- i

        if(is.null(hapmap36)){
            hapmap36 <- tmp
        }else{
            hapmap36 <- merge(hapmap36, tmp, all=TRUE)
        }
    }

}

rutgers36$Primer.SNP_ref_name <- NULL
names(rutgers36) <- c('marker', 'type', 'n.meioses', 'heterozygosity',
                      'phys.pos', 'cM.female', 'cM.male', 'cM', 'chr')
rutgers36 <- rutgers36[order(rutgers36$chr, rutgers36$phys.pos),]

rutgers37 <- subset(rutgers37, select = c(1:2, 4:9, 16))
names(rutgers37) <- c('marker', 'type', 'n.meioses', 'heterozygosity',
                      'phys.pos', 'cM', 'cM.female', 'cM.male', 'chr')
rutgers37 <- rutgers37[order(rutgers37$chr, rutgers37$phys.pos),]

names(hapmap36) <- c('phys.pos', 'cM', 'chr')
hapmap36 <- hapmap36[order(hapmap36$chr, hapmap36$phys.pos),]

save(rutgers36, file='../data/Maps/rutgers36.RData')
save(rutgers37, file='../data/Maps/rutgers37.RData')
save(hapmap36, file = '../data/Maps/hapmap36.RData')
