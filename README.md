ALDdata
=======

Data for use in analysis of Admixture Linkage Disequilibrium (see http://github.com/johnsonra/ALDdata for details).

This is an R package that is suggested for use with ALDsuite

File size is too big to upload complete the R package, so 
users will need to build their own local copy. Download a zipped 
archive of the repository, unzip and enter the following 
command in the terminal window from within the main directory.

sudo R CMD ISNTALL ALDdata


Versions
========

1.3.0
  - Added phased data for ASW

1.2.0
  - Added variant allele frequency for all SNPs in YRI, CEU, CHB and JPT
  - Added ALDdata to github

1.1.1
  - Fixed minor data bugs

1.1.0
  - Added CHB and JPT
  - Fixed cM measures on the X chromosome in rutgers36 and rutgers37 data sets

1.0.0
  - Initial split of data into it's own package