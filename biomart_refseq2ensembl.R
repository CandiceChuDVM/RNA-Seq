# If you are using UCSC genome and annotation with Tuxedo suite tools, 
# you might need to converse your refseq ID into Ensembl Gene ID fisrt
# then feed it into the Panther Classification System for gene ontology analysis.

# This R script provides useful solution to ID conversion.

# Install BioMart
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")

# Enter BioMart
library("biomaRt")

# Choose the database
ensembl=useMart("ensembl")

# Select canine genome in the Ensembl database
ensembl = useDataset("cfamiliaris_gene_ensembl",mart=ensembl)

# Import you input refseq.list into R. 
# The refseq.list can be made from editing your Cuffdiff output gene_exp.diff. 
# Use the following commands in linux to get the ideal refseq.list:
# $less gene_exp.diff | grep yes | grep rapid | grep control | cut -f3 | grep -v "-" > refseq.list.txt
# If you wish you have numeric values for Panther, please do change cut -f3 into cut -f3,x 
# (x means the order of column that you want to put into refseq.list)

mydata = read.table("refseq.list.txt") 

# Convert refseq to Ensembl Gene ID
results<- getBM(attributes = c("refseq_mrna","ensembl_gene_id"), filters="refseq_mrna", values=mydata, mart=ensembl)

# Export your Ensembl Gene ID list for Panther
write.table(results[,2], file="mydata.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Check mydata.txt in your working directory. Now you have a proper list ID for Panther!