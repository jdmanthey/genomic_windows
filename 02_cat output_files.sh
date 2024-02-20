# in windows directory

# combine the output for different analyses into a single file each
# first add a header for each file (repeat for each stat you want to keep)
grep 'pop1' CHR_10__10000001__10050000__stats.txt > ../window_heterozygosity.txt


# add the relevant stats to each file (repeat for each stat you want to keep)
for i in $( ls *txt ); do grep 'heterozygosity' $i >> ../window_heterozygosity.txt; done




# in R 
# in directory above windows directory
# to combine trees to single output file with another file about locations of each tree
options(scipen=999)

# list all the files in the trees directory
x_files <- list.files("windows", full.names=T)

# find the chromosome, start, and end for each tree
x_names <- list.files("windows")
x_chrom <- sapply(strsplit(sapply(strsplit(x_names, "RAxML_bipartitions."), "[[", 2), "__"), "[[", 1)
x_start <- sapply(strsplit(x_names, "__"), "[[", 2)
x_end <- sapply(strsplit(sapply(strsplit(x_names, "__"), "[[", 3), ".tre"), "[[", 1)

# write tree info
write.table(cbind(x_chrom, x_start, x_end), file="output_name_tree_info.txt", sep="\t", quote=F, row.names=F, col.names=F)

# trees into one file
tree_list <- list()
for(a in 1:length(x_files)) {
	tree_list[[a]] <- scan(x_files[a], what="character")
}
tree_list <- unlist(tree_list)
write(tree_list, file="output_name.trees", ncolumns=1)
