# hapmap.R
# Creation of hapmap datasets for use in mald package
# Randall Johnson
# Leidos Biomedical Research, Inc


library(ALDsuite)
library(ALDdata)
source('../R/helpers.R')

#############################
# Phased and LD data import #
#############################

# populations
ceu <- list()
yri <- list()
jpt <- list()
chb <- list()

asw <- list()

# for the final map
rs <- character()
chr <- numeric()
pos <- numeric()
ref <- character()
var <- character()
f.ceu <- numeric()
f.yri <- numeric()
f.jpt <- numeric()
f.chb <- numeric()

# ids for CHB and JPT
chb.ids <- sapply(strsplit(as.character(read.table('../data/CHB/pedinfo2sample_CHB.txt')$V6),
                           ":", fixed = TRUE), `[`, 5)
jpt.ids <- sapply(strsplit(as.character(read.table('../data/JPT/pedinfo2sample_JPT.txt')$V6),
                           ":", fixed = TRUE), `[`, 5)

for(i in 1:23)
{
    ### YRI ###
    phased <- get.phased(chr = i, pop = 'YRI')
    LD <- get.LD(chr = i, pop = 'YRI', keep = phased$rsID)

    # collect information
    rs.tmp <- phased$rsID
    chr.tmp <- rep(i, dim(phased)[1])
    pos.tmp <- phased$position_b36

    # get rid of unnecessary information
    phased$rsID <- NULL
    phased$position_b36 <- NULL
    phased$phys_position <- NULL
    rownames(phased) <- rs.tmp

    LD$fbin <- NULL

    tabs <- apply(phased, 1, table)

    # collect reference variables
    ref.tmp <- sapply(tabs, function(x) names(x)[order(-x)[1]])
    var.tmp <- sapply(tabs, function(x) names(x)[order(-x)[2]])

    drop.me <- numeric()

    # convert to 0/1 values
    for(j in 1:length(phased))
    {
        if(i == 23 & is.na(phased[[j]][1]))
        {
            drop.me <- c(drop.me, j)
        }else{
            phased[[j]] <- as.numeric(phased[[j]] != ref.tmp)
        }
    }

    if(i == 23 & length(drop.me) > 0)
        phased <- subset(phased, select = (1:length(phased))[-drop.me])
    phased <- t(phased)

    save(phased, LD, file = paste('../data/YRI/bin/yri', i, '.RData', sep = ''))
    f.yri <- c(f.yri, apply(phased, 2, mean))


    ### ASW ###
    if(i == 20)
    {
        phased <- get.phased(chr = i, pop = 'ASW')

        # get rid of unnecessary information
        phased$rsID <- NULL
        phased$position_b36 <- NULL
        phased$phys_position <- NULL
        rownames(phased) <- rs.tmp

        # find any variant alleles that don't exist in YRI (because there are no hets)
        to.check <- which(is.na(var.tmp))
        if(length(tabs) > 0)
        {
            tabs <- apply(phased[to.check,], 1, table)

            for(j in 1:length(to.check))
            {
                if(any(names(tabs[[j]]) != ref.tmp[to.check[j]]))
                    var.tmp[to.check[j]] <- names(tabs[[j]])[names(tabs[[j]]) != ref.tmp[to.check[j]]]
            }
        }

        drop.me <- numeric()

        # convert to 0/1 values
        for(j in 1:length(phased))
        {
            if(i == 23 & is.na(phased[[j]][1]))
            {
                drop.me <- c(drop.me, j)
            }else{
                phased[[j]] <- as.numeric(phased[[j]] != ref.tmp)
            }
        }

        if(i == 23 & length(drop.me) > 0)
            phased <- subset(phased, select = (1:length(phased))[-drop.me])
        phased <- t(phased)

        save(phased, file = paste('../Data/ASW/bin/asw', i, '.RData', sep = ''))
    }


    ### CEU ###
    phased <- get.phased(chr = i, pop = 'CEU')
    LD <- get.LD(chr = i, pop = 'CEU', keep = phased$rsID)

    # double check we are sorted correctly
    if(any(phased$rsID != rs.tmp))
        stop("Names mismatch on chromosome ", i)

    # get rid of unnecessary information
    phased$rsID <- NULL
    phased$position_b36 <- NULL
    phased$phys_position <- NULL
    rownames(phased) <- rs.tmp

    LD$fbin <- NULL

    # find any variant alleles that don't exist in YRI (because there are no hets)
    to.check <- which(is.na(var.tmp))
    if(length(tabs) > 0)
    {
        tabs <- apply(phased[to.check,], 1, table)

        for(j in 1:length(to.check))
        {
            if(any(names(tabs[[j]]) != ref.tmp[to.check[j]]))
                var.tmp[to.check[j]] <- names(tabs[[j]])[names(tabs[[j]]) != ref.tmp[to.check[j]]]
        }
    }

    drop.me <- numeric()

    # convert to 0/1 values
    for(j in 1:length(phased))
    {
        if(i == 23 & is.na(phased[[j]][1]))
        {
            drop.me <- c(drop.me, j)
        }else{
            phased[[j]] <- as.numeric(phased[[j]] != ref.tmp)
        }
    }

    if(i == 23 & length(drop.me) > 0)
        phased <- subset(phased, select = (1:length(phased))[-drop.me])
    phased <- t(phased)

    save(phased, LD, file = paste('../data/CEU/bin/ceu', i, '.RData', sep = ''))
    f.ceu <- c(f.ceu, apply(phased, 2, mean))


    ### CHB ###
    if(i == 23)
    {
        phased <- get.phased(chr = 23, pop = 'CHB')
    }else{
        phased <- get.phased(chr = i, pop = 'JPT+CHB')
        phased <- subset(phased, select = c('rsID', 'position_b36',
                                            names(phased)[substr(names(phased), 1, 7) %in% chb.ids]))
    }

    # order of names get messed up here :P
    rs.ord <- 1:length(rs.tmp)
    names(rs.ord) <- rs.tmp

    phased <- phased[order(rs.ord[phased$rsID]),]


    LD <- get.LD(chr = i, pop = 'CHB', keep = phased$rsID)

    # double check we are sorted correctly
    if(any(phased$rsID != rs.tmp))
        stop("Names mismatch on chromosome ", i)

    # get rid of unnecessary information
    phased$rsID <- NULL
    phased$position_b36 <- NULL
    phased$phys_position <- NULL
    rownames(phased) <- rs.tmp

    LD$fbin <- NULL

    # find any variant alleles that don't exist in YRI (because there are no hets)
    to.check <- which(is.na(var.tmp))
    if(length(tabs) > 0)
    {
        tabs <- apply(phased[to.check,], 1, table)

        for(j in 1:length(to.check))
        {
            if(any(names(tabs[[j]]) != ref.tmp[to.check[j]]))
                var.tmp[to.check[j]] <- names(tabs[[j]])[names(tabs[[j]]) != ref.tmp[to.check[j]]]
        }
    }

    drop.me <- numeric()

    # convert to 0/1 values
    for(j in 1:length(phased))
    {
        if(i == 23 & is.na(phased[[j]][1]))
        {
            drop.me <- c(drop.me, j)
        }else{
            phased[[j]] <- as.numeric(phased[[j]] != ref.tmp)
        }
    }

    if(i == 23 & length(drop.me) > 0)
        phased <- subset(phased, select = (1:length(phased))[-drop.me])
    phased <- t(phased)

    save(phased, LD, file = paste('../data/CHB/bin/chb', i, '.RData', sep = ''))
    f.chb <- c(f.chb, apply(phased, 2, mean))


    ### JPT ###
    if(i == 23)
    {
        phased <- get.phased(chr = 23, pop = 'JPT')
    }else{
        phased <- get.phased(chr = i, pop = 'JPT+CHB')
        phased <- subset(phased, select = c('rsID', 'position_b36',
                                            names(phased)[substr(names(phased), 1, 7) %in% jpt.ids]))
    }

    # order of names get messed up here :P
    rs.ord <- 1:length(rs.tmp)
    names(rs.ord) <- rs.tmp

    phased <- phased[order(rs.ord[phased$rsID]),]

    LD <- get.LD(chr = i, pop = 'JPT', keep = phased$rsID)

    # double check we are sorted correctly
    if(any(phased$rsID != rs.tmp))
        stop("Names mismatch on chromosome ", i)

    # get rid of unnecessary information
    phased$rsID <- NULL
    phased$position_b36 <- NULL
    phased$phys_position <- NULL
    rownames(phased) <- rs.tmp

    LD$fbin <- NULL

    # find any variant alleles that don't exist in YRI (because there are no hets)
    to.check <- which(is.na(var.tmp))
    if(length(tabs) > 0)
    {
        tabs <- apply(phased[to.check,], 1, table)

        for(j in 1:length(to.check))
        {
            if(any(names(tabs[[j]]) != ref.tmp[to.check[j]]))
                var.tmp[to.check[j]] <- names(tabs[[j]])[names(tabs[[j]]) != ref.tmp[to.check[j]]]
        }
    }

    drop.me <- numeric()

    # convert to 0/1 values
    for(j in 1:length(phased))
    {
        if(i == 23 & is.na(phased[[j]][1]))
        {
            drop.me <- c(drop.me, j)
        }else{
            phased[[j]] <- as.numeric(phased[[j]] != ref.tmp)
        }
    }

    if(i == 23 & length(drop.me) > 0)
        phased <- subset(phased, select = (1:length(phased))[-drop.me])
    phased <- t(phased)

    save(phased, LD, file = paste('../data/JPT/bin/jpt', i, '.RData', sep = ''))
    f.jpt <- c(f.jpt, apply(phased, 2, mean))


    rs <- c(rs, rs.tmp)
    chr <- c(chr, chr.tmp)
    pos <- c(pos, pos.tmp)
    ref <- c(ref, ref.tmp)
    var <- c(var, var.tmp)
}

save.image('tmp.RData')

#########################
# Creation of prior map #
#########################

hapmap <- data.frame(rs = rs,
                     chr = chr,
                     pos = pos,
                     ref = ref, # to be coded as 0
                     var = var, # to be coded as 1
                     f.yri = f.yri,
                     f.ceu = f.ceu,
                     f.chb = f.chb,
                     f.jpt = f.jpt,
                     stringsAsFactors = FALSE)

hapmap <- hapmap[order(hapmap$chr, hapmap$pos),]

hapmap$cM <- gen.calc(hapmap$chr, hapmap$pos, n.extrap = 20, map = "rutgers36")$gen.pos


save(hapmap, file = '../data/hapmap.RData')


##################################
# These didn't end up sorted! :P #
##################################

for(i in 1:23)
{
    hapmap.sub <- subset(hapmap, chr == i)

    for(pop in c('CEU', 'YRI', 'CHB', 'JPT'))
    {
        load(paste('../data/', pop, '/bin/', tolower(pop), i, '.RData', sep = ''))
        phased <- phased[,hapmap.sub$rs]

        save(phased, LD, file = paste('../data/', pop, '/bin/', tolower(pop), i, '.RData', sep = ''))
    }
}
