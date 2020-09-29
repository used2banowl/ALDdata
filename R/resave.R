# resave.R
# re-save .RData files to get better compression (run from root directory of repository)

files <- paste0('data/', system('ls data | grep .RData', intern = TRUE))

for(f in files)
{
  if(f %in% c('data/datalist', 'data/hapmap.RData', 
              'data/hapmap36.RData', 'data/rutgers36.RData',
              'data/rutgers37.RData'))
  {
    next()
  }else if(length(grep('asw', f)) == 1){
    next()
    load(f)
    save(phased, file = f, compress = 'xz')
    rm(phased)
  }else{
    load(f)
    save(phased, LD, file = f, compress = 'xz')
    rm(phased, LD)
  }
}

load('data/hapmap.RData')
save(hapmap, file = 'data/hapmap.RData', compress = 'xz')

load('data/hapmap36.RData')
save(hapmap36, file = 'data/hapmap36.RData', compress = 'xz')

load('data/rutgers36.RData')
save(rutgers36, file = 'data/rutgers36.RData', compress = 'xz')

load('data/rutgers37.RData')
save(rutgers37, file = 'data/rutgers37.RData', compress = 'xz')
