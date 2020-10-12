##### Determination of Assembly Processes #####
### By Stacey J. Doherty ###
### last modified 8 October 2020 ###

## This script was used to identify assembly processes structuring bacterial communities using 16S rRNA data.
## This script is intended to accompany the publication Doherty et al. 2020. "The transition from stochastic to deterministic bacterial community assembly during permafrost thaw succession"
## The original code was published by Stegen et al. 2013. "Quantifying community assembly processes and identifying features that impose them." ISME J.
## Please cite the original paper if using this script! I have merely added notes for clarification. Orginial code can be found on Stegen GitHub.


##### This section will generate the bNTI pairwise matrix for the dataset #####

## Change to the directory on your computer that contains the ASV table and associated phylogeny
## I used the rarefied ASV table and phylogeny. This one phylogeny will include all ASVs present in your dataset.
## Note that the 'slash' needs to be a forward slash like this /
setwd("/path_to_directory_containing_files")

## Load the picante library (Kembel et al., 2010)
## if not already installed, use install.packages('picante')
library(picante)

## Read in ASV table. It should be formatted as a .csv file. Rows are ASVs and columns are samples.
## This table was generated in qiime2, exported to a .biom file, then converted to a .csv file (code not shown).
## Note, Stegen et al. (2013) performed this analysis with OTU data, but for the sake of this code OTU = ASV.
otu = read.csv("feature-table.csv",header=T,row.names=1);
dim(otu); # this gives the dimensions
otu[1:5,1:5]; # this gives a look at the first 5 rows and columns

## Read in the phylogeny. It should be formatted as a .nwk file. 
## Phylogeny was constructed in qiime2, exported as a .nwk format
phylo = read.tree("tree.nwk");
phylo; # a summary of the phylogeny
plot.phylo(phylo,typ="fan"); # a quick plot, this may not work well depending on computer specs, but it isn't necessary to run.

## Make sure the names on the phylogeny are ordered the same as the names in ASV table
match.phylo.otu = match.phylo.data(phylo, otu);
str(match.phylo.otu);

## Calculate empirical betaMNTD

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
write.csv(beta.mntd.weighted,'betaMNTD_weighted.csv',quote=F);

identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

# Calculate randomized betaMNTD

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"weighted_bNTI.csv",quote=F);

pdf("weighted_bNTI_Histogram.pdf")
hist(weighted.bNTI)
dev.off()

## The weighted_bNTI.csv contains a pairwise matrix of bNTI values for the dataset.
## I opened the file in excel to assess and document assembly processes. 
## -2 < bNTI < 2 indicates stochastic process which can be further refined using values from RCbray in next section.
## bNTI < -2 indicates homogeneous selection
## bNTI > 2 indicates heterogeneous selection



##### This section will generate the RCbray pairwise matrix for the dataset #####

## This script is quite computationally heavy, I recommend running it on a super computer. 

## Change to the directory on your computer that contains the ASV table and associated phylogeny
## I used the rarefied ASV table and phylogeny. This one phylogeny will include all ASVs present in your dataset.
## Note that the 'slash' needs to be a forward slash like this /
setwd("/path_to_directory_containing_files")

## Read in ASV table 
## This table was generated in qiime2, exported to a .biom file, then converted to a .csv file (code not shown)
asv = read.csv("feature-table.csv",header=T,row.names=1);
asv = t(asv)
dim(asv); # this gives the dimensions
asv[1:5,1:5]; # this gives a look at the first 5 rows and columns

spXsite1 = asv

# Load vegan package 
library(vegan)

## Only deviation from Stegen code: SJD edited this function to use vegan for distance matrix creation.
raup_crick_abundance = function(spXsite, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE){
  
  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  
  ##By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). 
  ##Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. 
  ##The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number 
  ##of shared species to the calculation- this is highly recommended.  
  ## The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  
  ## Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  
  ## If ties are split (as we recommend) the dissimilarity (default) and similarity 
  ## (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) 
  ## or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). 
  ## If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. 
  ## The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  
  ## set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.  
  
  
  ## Note that the choice of how many plots (rows) to include has a real impact on the metric,
  ## as species and their occurrence frequencies across the set of plots is used to determine gamma and 
  ## the frequency with which each species is drawn from the null model	
  
  
  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }
  
  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  
  ##make the spXsite matrix into a new, pres/abs. matrix:
  ceiling(spXsite/max(spXsite))->spXsite.inc
  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite.inc, MARGIN=2, FUN=sum)
  
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(spXsite, MARGIN=2, FUN=sum)
  
  ##make_null:
  
  ##looping over each pairwise community combination:
  
  for(null.one in 1:(nrow(spXsite)-1)){
    for(null.two in (null.one+1):nrow(spXsite)){
      
      null_bray_curtis<-NULL
      for(i in 1:reps){
        
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(spXsite[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');			
        
        ##same for com2:
        com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(spXsite[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        
        null.spXsite = rbind(com1,com2); # null.spXsite;
        
        ##calculate null bray curtis
        null_bray_curtis[i] = vegdist(null.spXsite,method='bray');
        
      }; # end reps loop
      
      ## empirically observed bray curtis
      obs.bray = vegdist(spXsite[c(null.one,null.two),],method='bray');
      
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      
      rc = (num_less_than_in_null )/reps; # rc;
      
      if(split_ties){
        
        rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      };
      
      
      if(!classic_metric){
        
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        
        rc = (rc-.5)*2
      };
      
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      
      print(c(null.one,null.two,date()));
      
    }; ## end null.two loop
    
  }; ## end null.one loop
  
  if(as.distance.matrix){ ## return as distance matrix if so desired
    results<-as.dist(results)
  }	
  
  return(results)
  
}; ## end function

raup_crick_abundance(spXsite1, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE)

## The RCbray pairwise matrix will print to the terminal (console). This code could be adapted to write the results to a .csv file.
## I copied the results into a .csv file. The file was opened in excel to further refine the stochastic assembly processes. 
## -2 < bNTI < 2 and RCbray < -0.95 indicates homogenizing dispersal
## -2 < bNTI < 2 and RCbray > 0.95 indicates dispersal limitation and drift
## -2 < bNTI < 2 and -0.95 < RCbray < 0.95 indicates drift

