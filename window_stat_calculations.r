#functions used to calculate differentiation between two populations or diversity statistics within populations
# make sure the output file is already written in format (header):
# population1	population2	statistic	chromosome	start	end	number_sites	number_variable_sites	calculated_stat
# e.g., write(c("pop1", "pop2", "stat", "chr", "start", "end", "number_sites", "number_variable_sites", "calculated_stat"), ncolumns=9, file=outname, sep="\t")

# input is simplified vcf (entire vcf) and 3-column popmap (ordered the same as the vcf) and output file name
# and the input simple vcf file name that contains the chr and start and end information
# only calculate for invariant sites and biallelic snps
heterozygosity <- function(xxx, popmap, outname, filename) {
  # remove sites that are not either invariant or bi-allelic SNPs
  xxx <- xxx[nchar(xxx[,2]) == 1 & nchar(xxx[,3]) == 1, ]
  for(a in 1:nrow(popmap)) {
    output_rep <- c(popmap[a,1], "none", "heterozygosity", strsplit(filename, "__")[[1]][1], 
                    as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                    as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]))
    # select this individual
    a_rep <- xxx[,a+3]
    # remove missing genotypes
    a_rep <- a_rep[a_rep != "./."]
    # count number of sites
    a_total <- length(a_rep)
    # remove phasing information
    a_rep <- gsub("\\|", "/", a_rep)
    # find number of heterozygous sites
    a_het <- length(a_rep[a_rep == "0/1"])
    # add to output
    output_rep <- c(output_rep, a_total, a_het, a_het/a_total)
    write(output_rep, file=outname, append=T, ncolumns=9, sep="\t")
  }
}



# fis helper function to get expected heterozygosity for each snp in populations with 5 or more individuals
# to be applied across rows (i.e., each SNP)
# unbiased expected heterozygosity from Hahn book eqn. 3.1
expected_heterozygosity <- function(xxx) {
  # remove phasing information
  xxx <- gsub("\\|", "/", xxx)
  # check that there are at least 5 individuals
  if(length(xxx) > 4) {
    # remove missing genotypes
    xxx <- xxx[xxx != "./."]
    # sample size
    sample_size <- length(xx) * 2
    # allele freq
    p <- (length(xxx[xxx == "0/0"]) * 2 + length(xxx[xxx == "0/1"])) / sample_size
    # expected heterozygosity = 2pq (with sample size correction)
    eh <- sample_size / (sample_size - 1) * (1 - p^2 - (p-1)^2)
    # output
    return(eh)
  } else {
    return(NA)
  }
}



# fis calculation, same input as heterozygosity above 
# population requires 5 individuals to calculate
# Fis from Nei 1977 eqn. 14: F-statistics and analysis of gene diversity in subdivided populations 
# unbiased expected heterozygosity from Hahn book eqn. 3.1
fis <- function(xxx, popmap, outname, filename) {
  # remove sites that are not either invariant or bi-allelic SNPs
  xxx <- xxx[nchar(xxx[,2]) == 1 & nchar(xxx[,3]) == 1, ]
  for(a in 1:nrow(popmap)) {
    output_rep <- c(popmap[a,1], "none", "fis", strsplit(filename, "__")[[1]][1], 
                    as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                    as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]))
    # select this individual
    a_rep <- xxx[,a+3]
    # select all individuals in the same population
    a_pop <- xxx[,4:ncol(xxx)]
    a_pop <- a_pop[ ,popmap$PopName == popmap$PopName[a]]
    if(is.null(nrow(a_pop)) == FALSE) {
    	    # calculate expected heterozygosity for each site
    	a_eh <- apply(a_pop, 1, expected_heterozygosity)
    	# average expected heterozygosity across all sites
    	a_eh <- mean(a_eh)
    	# calculate observed heterozygosity for this individual
    	# remove missing genotypes
    	a_rep <- a_rep[a_rep != "./."]
    	# remove phasing information
    	a_rep <- gsub("\\|", "/", a_rep)
    	# find observed heterozygosity
    	a_het <- length(a_rep[a_rep == "0/1"]) / length(a_rep)
    	# compare expected and observed heterozygosities to calculate Fis
    	fis <- 1 - (a_eh - a_het) / a_eh
    	# add to output
    	output_rep <- c(output_rep, length(a_rep), length(a_rep[a_rep == "0/1"]), fis)
    	write(output_rep, file=outname, append=T, ncolumns=9, sep="\t")
    }
  }
}



# function to determine if a site is polymorphic (to be applied across rows)
# used in pi_tajima_theta function
polymorphic_function <- function(xxx) {
  return(length(unique(xxx)))
}



# function to determine frequency of allele 1 (to be applied across rows)
# used in pi_tajima_theta function
p_function <- function(xxx) {
  xxxx <- c(sapply(strsplit(xxx, "/"), "[[", 1), sapply(strsplit(xxx, "/"), "[[", 2))
  xxxx <- length(xxxx[xxxx == 0])
  return(xxxx)
}



# function to determine frequency of allele 2 (to be applied across rows)
# used in pi_tajima_theta function
q_function <- function(xxx) {
  xxxx <- c(sapply(strsplit(xxx, "/"), "[[", 1), sapply(strsplit(xxx, "/"), "[[", 2))
  xxxx <- length(xxxx[xxxx == 1])
  return(xxxx)
}



# input is simplified vcf (subsampled to single population (plus initial three columns)), 
# the name of the populations, and output file name
# and the input simple vcf file name that contains the chr and start and end information
# only calculate for invariant sites and biallelic snps
# no missing data allowed
# estimate watterson's theta, nucleotide diversity, and tajima's d with the formulae from 
# Carlson et al. 2005: 10.1101/gr.4326505 (equations 1, 2, and 4)
# e.g., xxx <- test[,c(1,2,3,5,6,7)]
pi_tajima_theta <- function(xxx, popname, outname, filename) {
  # remove sites that are not either invariant or bi-allelic SNPs
  xxx <- xxx[nchar(xxx[,2]) == 1 & nchar(xxx[,3]) == 1, ]
  # remove phasing information
  for(a in 4:ncol(xxx)) {
    xxx[,a] <- gsub("\\|", "/", xxx[,a])
  }
  # remove sites with missing data
  for(a in 4:ncol(xxx)) {
    xxx <- xxx[xxx[,a] != "./.",]
  }
  
  # define number of chromosomes sampled
  n_chromosomes <- (ncol(xxx) - 3) * 2
  # define number of sites genotyped for this population without missing data
  n_sites <- nrow(xxx)
  
  # continue if there are sites genotyped
  if(n_sites > 0) {
    # subset to only polymorphic sites
    xxx <- xxx[apply(xxx[,4:ncol(xxx)], 1, polymorphic_function) > 1, ]
    
    # count polymorphic sites
    n_polymorphic <- nrow(xxx)
    
    # determine denominator of theta calculation
    theta_denominator <- sum(1 / seq(from=1, to=(n_chromosomes - 1)))
    
    # continue if there are polymorphic sites
    if(n_polymorphic > 0) {
      # calculate per site theta
      theta <- (n_polymorphic / theta_denominator) / n_sites
      
      # calculate frequencies of alleles 1 and 2 for each polymorphic site
      p_freq <- as.vector(apply(xxx[,4:ncol(xxx)], 1, p_function)) / n_chromosomes
      q_freq <- as.vector(apply(xxx[,4:ncol(xxx)], 1, q_function)) / n_chromosomes
      
      # calculate per site nucleotide diversity (pi)
      pi <- sum(2 * p_freq * q_freq) * (n_chromosomes / (n_chromosomes - 1)) / n_sites
      
      # calculate D using the variance of d (and all necessary subcomponents)
      a1 <- sum(1 / seq(from=1, to=(n_chromosomes - 1)))
      a2 <- sum(1 / (seq(from=1, to=(n_chromosomes - 1))^2))
      b1 <- (n_chromosomes + 1) / (3 * (n_chromosomes - 1))
      b2 <- (2 * (n_chromosomes^2 + n_chromosomes + 3)) / (9 * n_chromosomes * (n_chromosomes - 1))
      c1 <- b1 - (1 / a1)
      c2 <- b2 - ((n_chromosomes + 2) / (a1 * n_chromosomes)) + (a2 / (a1^2))
      e1 <- c1 / a1
      e2 <- c2 / ((a1^2) + a2)
      Var_d <- (e1 * n_polymorphic) + (e2 * n_polymorphic * (n_polymorphic - 1))
      Tajima_D <- ((pi * n_sites - theta * n_sites)) / sqrt(Var_d)
      
      # write output
      output_rep1 <- c(popname, "none", "theta", strsplit(filename, "__")[[1]][1], 
                       as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                       as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]),
                       n_sites, n_polymorphic, theta)
      output_rep2 <- c(popname, "none", "pi", strsplit(filename, "__")[[1]][1], 
                       as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                       as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]),
                       n_sites, n_polymorphic, pi)
      output_rep3 <- c(popname, "none", "Tajima_D", strsplit(filename, "__")[[1]][1], 
                       as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                       as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]),
                       n_sites, n_polymorphic, Tajima_D)
      write(output_rep1, file=outname, append=T, ncolumns=9, sep="\t")
      write(output_rep2, file=outname, append=T, ncolumns=9, sep="\t")
      write(output_rep3, file=outname, append=T, ncolumns=9, sep="\t")
    } else { # if no polymorphic sites retained
      output_rep1 <- output_rep1 <- c(popname, "none", "theta", strsplit(filename, "__")[[1]][1], 
                                      as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                                      as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]),
                                      n_sites, n_polymorphic, "NA")
      output_rep2 <- output_rep1 <- c(popname, "none", "pi", strsplit(filename, "__")[[1]][1], 
                                      as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                                      as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]),
                                      n_sites, n_polymorphic, "NA")
      output_rep3 <- output_rep1 <- c(popname, "none", "Tajima_D", strsplit(filename, "__")[[1]][1], 
                                      as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                                      as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]),
                                      n_sites, n_polymorphic, "NA")
      write(output_rep1, file=outname, append=T, ncolumns=9, sep="\t")
      write(output_rep2, file=outname, append=T, ncolumns=9, sep="\t")
      write(output_rep3, file=outname, append=T, ncolumns=9, sep="\t")
    }
  } else { # if no sites retained
    output_rep1 <- output_rep1 <- c(popname, "none", "theta", strsplit(filename, "__")[[1]][1], 
                                    as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                                    as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]),
                                    n_sites, "NA", "NA")
    output_rep2 <- output_rep1 <- c(popname, "none", "pi", strsplit(filename, "__")[[1]][1], 
                                    as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                                    as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]),
                                    n_sites, "NA", "NA")
    output_rep3 <- output_rep1 <- c(popname, "none", "Tajima_D", strsplit(filename, "__")[[1]][1], 
                                    as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                                    as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]),
                                    n_sites, "NA", "NA")
    write(output_rep1, file=outname, append=T, ncolumns=9, sep="\t")
    write(output_rep2, file=outname, append=T, ncolumns=9, sep="\t")
    write(output_rep3, file=outname, append=T, ncolumns=9, sep="\t")
  }
  
}



# function to determine if a site is only missing data for a population (used across rows of vcf)
# used in differentiation function
total_missing <- function(xxx) {
  return(length(xxx[xxx != "./."]) > 0)
}	



# function to determine if a site is polymorphic (to be applied across rows) after removing missing
# used in differentiation function
polymorphic_function2 <- function(xxx) {
  xxx <- xxx[xxx != "./."]
  return(length(unique(xxx)))
}			



# input is two simplified vcfs (subsampled to single population), already filtered for invariant/biallelic SNPs, 
# the names of the populations, and output file name
# and the input simple vcf file name that contains the chr and start and end information
# only calculate for invariant sites and biallelic snps
# fst is calculation of Reich et al. 2009
# fst is the reich et al. 2009 estimator for small sample sizes
# equation presented nicer in Willing et al. 2012 page 9
# Dxy calculation is from Tavares et al. 2018: pnas.1801832115
differentiation <- function(xxx1, xxx2, popname1, popname2, outname, filename) {
  # remove phasing information
  for(a in 1:ncol(xxx1)) {
    xxx1[,a] <- gsub("\\|", "/", xxx1[,a])
  }
  for(a in 1:ncol(xxx2)) {
    xxx2[,a] <- gsub("\\|", "/", xxx2[,a])
  }
  
  # remove sites that are completely missing from either population
  keep1 <- apply(xxx1, 1, total_missing)
  keep2 <- apply(xxx2, 1, total_missing)
  xxx1 <- xxx1[keep1 == TRUE & keep2 == TRUE, ]
  xxx2 <- xxx2[keep1 == TRUE & keep2 == TRUE, ]
  
  # count the total number of included genotyped sites at this point
  n_sites <- nrow(xxx1)
  
  # combine the two matrices to find sites that are variant w/in and between the two pops
  xxx_combined <- cbind(xxx1, xxx2)
  variant_sites <- apply(xxx_combined, 1, polymorphic_function2)
  
  # keep only variant sites
  xxx1_variant <- xxx1[variant_sites > 1, ]
  xxx2_variant <- xxx2[variant_sites > 1, ]
  
  # count the number of variant sites
  n_variant_sites <- nrow(xxx1_variant)
  
  # loop for each polymorphic site to calculate dxy
  dxy_all <- list()
  for(a in 1:nrow(xxx1_variant)) {
    a_rep1 <- as.character(xxx1_variant[a,])
    a_rep2 <- as.character(xxx2_variant[a,])
    
    # remove missing
    a_rep1 <- a_rep1[a_rep1 != "./."]
    a_rep2 <- a_rep2[a_rep2 != "./."]
    
    # measure proportion of reference allele 
    a_ref1 <- (length(a_rep1[a_rep1 == "0/0"]) * 2 + length(a_rep1[a_rep1 == "0/1"]) * 1) / (length(a_rep1) * 2)
    a_ref2 <- (length(a_rep2[a_rep2 == "0/0"]) * 2 + length(a_rep2[a_rep2 == "0/1"]) * 1) / (length(a_rep2) * 2)
    
    # calc dxy
    dxy_all[[a]] <- a_ref1 * (1 - a_ref2) + a_ref2 * (1 - a_ref1)
  }
  dxy_all <- sum(unlist(dxy_all)) / n_sites
  
  
  # loop for each polymorphic site to calculate fst
  numerator_all <- list()
  denominator_all <- list()
  for(a in 1:nrow(xxx1_variant)) {
    a_rep1 <- as.character(xxx1_variant[a,])
    a_rep2 <- as.character(xxx2_variant[a,])
    
    # remove missing
    a_rep1 <- a_rep1[a_rep1 != "./."]
    a_rep2 <- a_rep2[a_rep2 != "./."]
    
    # number of individuals per population
    pop1_ind_count <- length(a_rep1) 
    pop2_ind_count <- length(a_rep2)
    
    # non-reference allele counts
    alt_allele_count1 <- (2 * length(a_rep1[a_rep1 == "1/1"]) + 1 * length(a_rep1[a_rep1 == "0/1"]))
    alt_allele_count2 <- (2 * length(a_rep2[a_rep2 == "1/1"]) + 1 * length(a_rep2[a_rep2 == "0/1"]))
    
    # total allele counts
    all_allele_count1 <- 2 * length(a_rep1)
    all_allele_count2 <- 2 * length(a_rep2)
    
    # expected heterozygosity for each population
    expected_het1 <- (alt_allele_count1 * (all_allele_count1 - alt_allele_count1)) / 
      (all_allele_count1 * (all_allele_count1 - 1))
    expected_het2 <- (alt_allele_count2 * (all_allele_count2 - alt_allele_count2)) / 
      (all_allele_count2 * (all_allele_count2 - 1))
    
    # find the fst numerator and denominator values for this snp (they all get summed and divided for 
    # the final estimate)
    numerator_all[[a]] <- (alt_allele_count1 / (2 * pop1_ind_count) - 
                             alt_allele_count2 / (2 * pop2_ind_count))^2 - (expected_het1 / (2 * pop1_ind_count)) - 
      (expected_het2 / (2 * pop2_ind_count))
    denominator_all[[a]] <- numerator_all[[a]] + expected_het1 + expected_het2		
  }
  # calculate total fst for this window
  fst_all <- sum(unlist(numerator_all)) / sum(unlist(denominator_all))
  
  # write to output for dxy and fst
  output_rep1 <- c(popname1, popname2, "Dxy", strsplit(filename, "__")[[1]][1], 
                   as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                   as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]),
                   n_sites, n_variant_sites, dxy_all)
  output_rep2 <- c(popname1, popname2, "Fst", strsplit(filename, "__")[[1]][1], 
                   as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                   as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]),
                   n_sites, n_variant_sites, fst_all)
  write(output_rep1, file=outname, append=T, ncolumns=9, sep="\t")
  write(output_rep2, file=outname, append=T, ncolumns=9, sep="\t")
  
}



# helper function to identify transition variants from biallelic alleles in rows (used in titv function)
transitions <- function(x1) {
  a_rep <- sort(x1)
  if(a_rep[1] == "A" & a_rep[2] == "G") {
    return(1)
  } else if(a_rep[1] == "C" & a_rep[2] == "T") {
    return(1)
  } else {
    return(0)
  }
}



# input is simplified vcf (entire vcf) and other files as in heterozygosity calcs, to calculate transition transversion ratio
titv <- function(xxx, popmap, outname, filename) {
  options(scipen=999)
  xxx <- xxx[nchar(xxx[,2]) == 1 & nchar(xxx[,3]) == 1, ]
  # filter to variable sites
  xxx <- xxx[xxx[,3] != ".", ]
  
  # sites used for titv 
  # total sites used same as het sites
  a_total <- nrow(xxx)
  a_het <- a_total
  
  # find number of transitions using helper function
  total_transitions <- sum(apply(xxx[,2:3], 1, transitions))
  # calculate transition / transversion ratio
  titv_ratio <- total_transitions / (a_total - total_transitions)
  
  # output info
  output_rep <- c("all_inds", "none", "titv", strsplit(filename, "__")[[1]][1], 
                   as.numeric(strsplit(strsplit(filename, "__")[[1]][2], "__")[[1]][1]),
                   as.numeric(strsplit(strsplit(filename, "__")[[1]][3], ".simple")[[1]][1]))
  output_rep <- output_rep <- c(output_rep, a_total, a_het, titv_ratio)
  write(output_rep, file=outname, append=T, ncolumns=9, sep="\t")
  
}



# input is one simplified vcfs, the popmap, and output file name
# and the input simple vcf file name that contains the chr and start and end information
# only for invariant sites and biallelic snps
# lastly, input number of sites needed to write file
create_fasta_from_vcf <- function(xxx, popmap, outname, filename, num_sites_needed) {
  options(scipen=999)
  # remove sites that are not either invariant or bi-allelic SNPs
  xxx <- xxx[nchar(xxx[,2]) == 1 & nchar(xxx[,3]) == 1, ]
  
  # keep going only if enough sites
  if(num_sites_needed <= nrow(xxx)) {
    # remove phasing information
    for(a in 4:ncol(xxx)) {
      xxx[,a] <- gsub("\\|", "/", xxx[,a])
    }
    
    # define names of individuals in output fasta
    output_names_fasta <- paste(">", popmap[,1], sep="")
    
    # subset the genotypes from the allele info
    allele_info <- xxx[,2:3]
    genotypes <- xxx[,4:ncol(xxx)]
    
    # convert all numbers in genotypes to actual bases and ambiguities
    for(a in 1:nrow(genotypes)) {
      if(allele_info[a,2] == ".") { # if non-polymorphic
        genotypes[a,] <- gsub("0/0", allele_info[a,1], genotypes[a,])
        genotypes[a,] <- gsub("\\./\\.", "?", genotypes[a,])
      } else { # if polymorphic
        both_genotypes <- sort(as.character(allele_info[a,]))
        if(both_genotypes[1] == "A" & both_genotypes[2] == "C") { het = c("A", "C") }
        if(both_genotypes[1] == "A" & both_genotypes[2] == "G") { het = c("A", "G") }
        if(both_genotypes[1] == "A" & both_genotypes[2] == "T") { het = c("A", "T") }
        if(both_genotypes[1] == "C" & both_genotypes[2] == "G") { het = c("C", "G") }
        if(both_genotypes[1] == "C" & both_genotypes[2] == "T") { het = c("C", "T") }
        if(both_genotypes[1] == "G" & both_genotypes[2] == "T") { het = c("G", "T") }
        genotypes[a,] <- gsub("0/0", allele_info[a,1], genotypes[a,])
        genotypes[a,] <- gsub("\\./\\.", "?", genotypes[a,])
        genotypes[a,] <- gsub("1/1", allele_info[a,2], genotypes[a,])
        genotypes[a,][genotypes[a,] == "0/1"] <- sample(het, length(genotypes[a,][genotypes[a,] == "0/1"]), replace=T)
      }
    }
    
    # write output
    for(a in 1:ncol(genotypes)) {
      if(a == 1) {
        write(output_names_fasta[a], file=outname, ncolumns=1)
        write(paste(genotypes[,a], collapse=""), file=outname, ncolumns=1, append=T)
      } else {
        write(output_names_fasta[a], file=outname, ncolumns=1, append=T)
        write(paste(genotypes[,a], collapse=""), file=outname, ncolumns=1, append=T)
      }
    }
  }
}




