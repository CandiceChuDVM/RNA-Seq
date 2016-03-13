All blockquotes are adapted from [chapkinlab](https://github.com/chapkinlab/sequencing-pipeline/wiki). For example:   

>This is a guide to get you started quickly with beginning to use our sequencing pipeline for analyzing your RNA-seq datasets. 

DMSO is the control group and TCDD is the treatment group. We will be using data related to only these groups. There are six samples in total. Ignore the rest. Sequencing data is from two lanes, and these have to be concatenated before processing.

#1. Installation and Setup  
##20160310 Get the files
Log in VPN  
Open terminal on mac  
Type `ssh -p 22 candice@nfsc-oracle.tamu.edu`  
Type password  
>All the pipeline programs are on Github. These files might need some editing to accommodate project-specific configuration of the pipeline. If you are setting up your pipeline new, navigate to your desired installation directory, and issue the following command.    

Under folder `Assignment`, exicute `git clone https://github.com/chapkinlab/sequencing-pipeline.git`    

>The pipeline makes use of several external open source tools that you'll need to have installed on your computer. These include:   
+The STAR RNA-Seq aligner  
+FastQC  
+Python (with the Pandas library)      

>You'll also need a reference genome and annotation. You can obtain that from [Ensembl](http://useast.ensembl.org/info/data/ftp/index.html).  

    $ wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz       
    $ wget ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz   
##20160311
>To verify that the pipeline is installed and up-to-date, navigate to the `sequencing-pipeline` directory, and issue a `git pull`. An `Already up-to-date` message confirms the update.  

#2. Obtaining Sequencing Datasets  
    $ `cp -r /mnt/nas/Organoid-data-04242014/FastqFile ./home/candice/Assignment/FastqFile`  

#3. Quality Check  
>If the sequencing facility has not done a quality check on the sequencing data, you can run `preqc.sh` to check on the quality of the sequencing.  

##20160312
    $ `mkdir FastqFiles_fastq`    
`gunzip` every fastq.gz files    
Create a experiment list in `Assignment/sequencing-pipeline/lists`: 
    $ `vim candice_list`    

In `candice_list`:  

    species="grch38-human"  
    samplelist=(\  
    /home/candice/Assignment/FastqFile/Sample_DMSO1/DMSO1_AGTTCC_L001_R1_001.fastq.gz \  
    /home/candice/Assignment/FastqFile/Sample_DMSO1/DMSO1_AGTTCC_L002_R1_001.fastq.gz \  
    /home/candice/Assignment/FastqFile/Sample_DMSO2/DMSO2_GTGAAA_L001_R1_001.fastq.gz \  
    /home/candice/Assignment/FastqFile/Sample_DMSO2/DMSO2_GTGAAA_L002_R1_001.fastq.gz \  
    /home/candice/Assignment/FastqFile/Sample_DMSO3/DMSO3_TGACCA_L001_R1_001.fastq.gz \  
    /home/candice/Assignment/FastqFile/Sample_DMSO3/DMSO3_TGACCA_L002_R1_001.fastq.gz \  
    /home/candice/Assignment/FastqFile/Sample_TCDD1/TCDD1_CCGTCC_L001_R1_001.fastq.gz \  
    /home/candice/Assignment/FastqFile/Sample_TCDD1/TCDD1_CCGTCC_L002_R1_001.fastq.gz \  
    /home/candice/Assignment/FastqFile/Sample_TCDD2/TCDD2_CGATGT_L001_R1_001.fastq.gz \  
    /home/candice/Assignment/FastqFile/Sample_TCDD2/TCDD2_CGATGT_L002_R1_001.fastq.gz \  
    /home/candice/Assignment/FastqFile/Sample_TCDD3/TCDD3_CAGATC_L001_R1_001.fastq.gz \  
    /home/candice/Assignment/FastqFile/Sample_TCDD3/TCDD3_CAGATC_L002_R1_001.fastq.gz \
    )    

From the sequencing pipeline directory, run `main-scripts/preqc.sh lists/candice_list`  

#4. Dataset Pre processing and List files  
>If your sequences were run on multiple lanes, you will first need to concatenate your files so each sample has only one file. To easily do this, use the `main-scripts/concat.sh` script available in the sequencing-pipeline repository. 
 
>The `concat.sh` script automatically finds separate `fastq.gz` files from a single experiment (split across several sequencing lanes) and then combines them together into one file that you can then use in the mapping step.   

    $ `cd sequencing-pipeline/main-scripts`  
    $ `vi concat.sh`   

>Modifying the cut commmand on [this line](https://github.com/chapkinlab/sequencing-pipeline/blob/80061158cfebb1dda2c7806779b53466573cc337/main-scripts/concat.sh#L34) from 1-4 to 1-5 in `concat.sh`.   

>Remove the # commenting character) [this line](https://github.com/chapkinlab/sequencing-pipeline/blob/3155b43c4877023fc7a6b5699b77dd42d9bc2389/main-scripts/concat.sh#L50) in your local copy of the concat.sh script and then run it.  

    $ ``mkdir FastqFile_processed``  
    $ `./sequencing-pipeline/main-scripts/concat.sh ./FastqFile/ ./FastqFile_processed`   

(sample list is optional)  

#5. Mapping reads to the genome  
>Now that you know that the sequences have passed QC, it is time to map them against a reference genome. In order to make the processing of our datafiles easier, we need to make an "experiment list" file which describes the location of the samples to be operated on, and the reference genome that they should be mapped against. You will need to make an experiment list file `candice_list` for your samples. 

The genome that you will reference against (e.g. `grch38-human`) should be the first line in your `candice_list`. `map.sh` is the script to handle mapping your samples against the reference genome. This script is also usually run from the sequencing pipeline directory.   
##20160313
Edit `candice_list`:  

    species="grch38-human"  
    samplelist=(\  
    /home/candice/Assignment/FastqFile_processed/DMSO1.fastq.bz2 \  
    /home/candice/Assignment/FastqFile_processed/DMSO2.fastq.bz2 \  
    /home/candice/Assignment/FastqFile_processed/DMSO3.fastq.bz2 \  
    /home/candice/Assignment/FastqFile_processed/TCDD1.fastq.bz2 \  
    /home/candice/Assignment/FastqFile_processed/TCDD2.fastq.bz2 \  
    /home/candice/Assignment/FastqFile_processed/TCDD3.fastq.bz2 \  
    )    

Run `map.sh`:   

     $ mv Homo_sapiens.GRCh38.84.gtf.gz FastqFile_processed/  
     $ mv Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz FastqFile_processed/  
     $ gunzip Homo_sapiens.GRCh38.84.gtf.gz   
     $ gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz   
     $ cd ../sequencing-pipeline/
     $ main-scripts/map.sh lists/candice_list 2> err.log | tee out.log  

>If there is a problem running the post ribosomal reads, julia needs to be set up in the home directory. If you are on nfsc-oracle do the following:
cd ~  
julia  
julia> Pkg.update()  
julia> Pkg.add("HDF5")  
julia> Pkg.add("JLD")   
This will setup the ~/.julia directory, update the julia main packages, and add the HDF5 and JLD packages necessary to run misc-scripts/ribo.jl, which does the post ribosomal processing.   

#6. Analysis Summarization  
>In order to tell how the mapping went, we would like to see an overview (summary). Often times, an experiment will have additional 'metadata' describing the experimental conditions of each sample. It is typically useful to have this information included alongside the sequencing summary information to determine if there are any patterns (perhaps all treatment samples have lower number of mapped reads). This is often difficult using only sample names.  

>To do this, we need to create a `key file` as a comma delimited file (.csv). You have to have the first column be your sample name with a title of `sample`. The other columns may indicate the different treatment groups, metadata etc. It is usually easiest to do this in Excel (use one of the existing keyfiles for guidance) and then saving the output as a .csv file. This key file will also be used when running edgeR for running statistical comparisons (differential expression) between groups.    

Create a key file:

    $ vim candice-key.csv

candice-key.csv:

     sample, type
     DMSO1, control
     DMSO2, control
     DMSO3, control
     TCDD1, treatment
     TCDD2, treatment
     TCDD3, treatment

    $ main-scripts/summary.py lists/candice_list keys/candice-key.csv

>This will create output in the analysis/<name of experiment list file> directory. In it you can find the raw counts of the experimental samples in both samples x genes and genes x samples format (in the <explist>-count.csv and <explist>-count.T.csv files respectively). You can also find a summary spreadsheet in <explist>.summary.csv. There is also an <elist>.h5 file which can be read efficiently if you are using Python using the Pandas library.   

>The columns of the `summary.csv` spreadsheet are as follows:
total-reads: Number of reads in the fastQ file   
uniq-reads: Number of unique reads in the FastQ file.    
grch38-reads: Number of mapped reads (counting each multi-mapper separately)   
grch38-uniq: Number of reads which mapped uniquely (once) to the reference genome (this is separate from uniq-reads above).   
grch38-multi: Number of reads which mapped multiple locations to the reference.   
grch38-annotated-reads: Number of resulting reads mapped to annotated regions of the genome (genes, lncrna, etc.)   
mito-reads: Number of reads mapping to mitochondrial regions of the reference.   
ercc-reads: Number of reads which mapped to ERCC transcript sequences.  
rrna-reads: " " ribosomal sequences.  
htseq-0-genes: Number of genes which had at least 1 read  
htseq-3-genes: Number of genes which had at least 4 reads  
htseq-10-genes: Number of genes which had at least 11 reads  

#7. Gene Differential Expression Analysis  
>Now you are ready to find the differentially expressed genes. You will need to edit `run-edger.R` to set up the comparisons you want to run.   
Edit `/home/candice/Assignment/sequencing-pipeline/main-scripts/run-edger.R`    
     $ vim run-edger.R   
press `esc` then type `:set nu` to display line numbers.    
>The list file and the key file locations need to be specified inside the script, in line 19. Replace this line with the one containing the actual key and list files to be used.    

      args = c("lists/candice_list","keys/candice-key.csv")   

>Then, provide an experiment name (result directory name) that is specific to the analysis, by modifying line 45.    

       ename = "edger-pair-treatment"    

>Add the keys (columns in the key file) across which the differential expression analysis has to performed. This is done by modifying line 46.    

        factors = key[order(rownames(key)), c("type")]     

>Here, "type" and "treatment" are the columns used. Replace these with the relevant ones in your key file. Also change line 49 with the same columns.    

        design = model.matrix(~type, data=factors)     

>The factor for the pairwise analysis is specified in line 50. "treatment" in this line has to be replaced with the factor against which the differential expression analysis is performed.    

        groups = factors$type  

>Run run-edger.R.    

      $ Rscript main-scripts/run_edger.R    

