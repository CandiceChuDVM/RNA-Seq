#
DMSO is the control group and TCDD is the treatment group. We will be using data related to only these groups. There are six samples in total. Ignore the rest. Sequencing data is from two lanes, and these have to be concatenated before processing.
#1. Installation and Setup  
##20160310 Get the files
Log in VPN  
Open terminal on mac  
Type `ssh -p 22 candice@nfsc-oracle.tamu.edu`  
Type password  
Under folder `Assignment`, exicute `git clone https://github.com/chapkinlab/sequencing-pipeline.git`  
`wget ftp://ftp.ensembl.org/pub/release-79/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz`  
`wget ftp://ftp.ensembl.org/pub/release-81/gtf/homo_sapiens/Homo_sapiens.GRCh38.81.gtf.gz`     
##20160311
To verify that the pipeline is installed and up-to-date, navigate to the `sequencing-pipeline` directory, and issue a `git pull`. An `Already up-to-date` message confirms the update.  
#2. Obtaining Sequencing Datasets  
#3. Quality Check  
#4. Dataset Pre processing and List files  
`cd sequencing-pipeline/main-scripts`  
`vi concat.sh`
modifying the cut commmand on [this line](https://github.com/chapkinlab/sequencing-pipeline/blob/80061158cfebb1dda2c7806779b53466573cc337/main-scripts/concat.sh#L34) from 1-4 to 1-5 in `concat.sh`.   
remove the # commenting character) [this line](https://github.com/chapkinlab/sequencing-pipeline/blob/3155b43c4877023fc7a6b5699b77dd42d9bc2389/main-scripts/concat.sh#L50) in your local copy of the concat.sh script and then run it.  
`mkdir FastqFile_processed`  
`./sequencing-pipeline/main-scripts/concat.sh ./FastqFile/ ./FastqFile_processed concat.sh_processed_list.txt`   
FastqFile/Sample_DMSO1/DMSO1_AGTTCC_L001_R1_001.fastq.gz
FastqFile/Sample_DMSO1/DMSO1_AGTTCC_L002_R1_002.fastq.gz
