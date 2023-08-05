---
title: "ChIP-seq2"
author: Shilpa sharma, Sili vettiyara sunil
date: 22.08.2022
theme: united
highlight: tango
output:
  html_document:
    toc: true
    toc_depth: 3
    number_sections: true
    toc_float:
      collapsed: true
      smooth_scroll: false
tags: chipseq
---


# The environment - task sheet
## 2 Get to know the system
1. Run sinfo -lN on the master node.
``` {bash}
sinfo -lN
```
How many nodes do you see? - 4
How many CPUs? - 12 
How much RAM does each node have? - 121000

![](https://i.imgur.com/K3KSP91.png)

2. Run df -h /vol/COMPEPIWS/. 
``` {bash}
df -h /vol/COMPEPIWS/.
```
What is the volume’s size? - 3.9T
![](https://i.imgur.com/0pcn2Kg.png)


3. Check the temp folder location by
``` {bash}
echo $TMPDIR: /vol/COMPEPIWS/tmp
```
4. Your group working directory is 
```
/vol/COMPEPIWS/groups/chipseq2
```
6. The raw data (fastq files) that you are going to use to run the pipelines is available under this folder:
```
/vol/COMPEPIWS/data/reduced/<assay>
```
![](https://i.imgur.com/uv6MHkW.png)





## 3 Tips and Tricks
``` {bash}
1. $ ssh -X -A -i OpenSSH_private_key [ELIXIR_LOGIN]@193.175.249.28
```    

``` {bash}
2. $ ssh -X -A  -i OpenSSH_private_key bibigrid-worker-1-2-etqcgckjrdsvsxe
```    

``` {bash}
3. $ screen -S screen_1
```

``` {bash}
4. Ctrl-a c
   Ctrl-a c
```
    
``` {bash}
5. Ctrl-a n
```
    
``` {bash}
6. Ctrl-a d
```
    
``` {bash}
7. Ctrl-d
```
    
``` {bash}
8. Ctrl-d
```
    
``` {bash}
9. $ ssh -X -A -i OpenSSH_private_key [ELIXIR_LOGIN]@193.175.249.28
```
    
``` {bash}
10. $ ssh -X -A  -i OpenSSH_private_key bibigrid-worker-1-2-etqcgckjrdsvsxe 
```

``` {bash}    
11. $ screen -ls
```
    
``` {bash}     
12. $ screen -r
```    

``` {bash}    
13. Ctrl-a n
``` 


    
## 4 Working with tables on the command line, or awk

3. Number of lines in this file: 671462

 ```{bash}
wc -l /vol/COMPEPIWS/pipelines/references/genome_genes.gtf
```  
    

4. Number of "exon" in this file (hint: grep, awk): 324748
 ``` {bash}
grep -o -i exon /vol/COMPEPIWS/pipelines/refrences/genome_genes.gtf | wc -l
```   

5. Number of "exon" that are longer than 1000bp considering cloumns 4 and 5 as coordinates of the exon :19593 
 ``` {bash}
awk '$3 = "exon" && $5-$4 > 1000' /vol/COMPEPIWS/pipelines/refrences/genome_genes.gtf > exon_genes.gtf
```     
``` {bash}
wc -l exon_genes.gtf 19593 exon_genes.gtf
```    
![](https://i.imgur.com/QVdneN2.png)

6. Number of "exon" that has Sox17 as gene_id: 20
```{bash}
grep -o -i Sox17 exon_genes.gtf | wc -l
```
    
![](https://i.imgur.com/2JZ64Gz.png)

7. Number of "exon" in chr2: 29373
``` {bash}
awk '$1 = "chr2"' /vol/COMPEPIWS/pipelines/references/genome_genes.gtf > chr2_exon.gtf
```

``` {bash}
wc -l chr2_exon.gtf
60606 chr2_exon.gtf
```
    
 ``` {bash}
awk '$3 = "exon"' chr2_exon.gtf | wc -l
```   
    
![](https://i.imgur.com/jUYWLWL.png)

8. Frequency of each feature in column 3: 286783 CDS, 324748 exon, 29978 start_codon, 29953 stop_codon
``` {bash}
awk -F '''{print $3}' /vol/COMPEPIWS/pipelines/references/genome_genes.gtf |sort | uniq -c
```
    
![](https://i.imgur.com/3YQZFl3.png)

9. Sort the file according to column1 and then column3. redirect the output into a new file
 ``` {bash}
sort -k 1 /vol/COMPEPIWS/pipelines/references/genome_genes.gtf > col1.gtf
```  
``` {bash}
sort -k 3 col1.gtf > col13.gtf
```    
    
``` {bash}
head col12.gtf
```
``` {bash}
head col13.gtf
```
![](https://i.imgur.com/q9SlliS.png)


## 5 Conda and Bioconda
a) Load the “base” environment if it is not loaded    
    
   ``` {bash}
source /vol/COMPEPIWS/conda/miniconda3/bin/activate

```
b) create a new conda environment under your group conda directory called as your group’s name
``` {bash}
conda create -p /vol/COMPEPIWS/groups/chipseq2/conda/chipseq2

``` 
c) activate the created environment
``` {bash}
conda activate /vol/COPMPEPIWS/groups/chipseq2/conda/chipseq2

``` 
d) install ggplot2
``` {bash}
conda install -c conda-forge r-ggplot2

``` 
install fastqc
``` {bash}
conda install -c bioconda fastqc

``` 
install bedtools
``` {bash}
conda install -c bioconda bedtools=2.22

``` 
    
e) make sure the tools have been installed properly by :
``` {bash}
$ fastqc --help
```
    
``` {bash}
$ bedtools --help
```
    
``` {bash}
$ R
```

f) update bedtools to the latest version, check the new version number
``` {bash}
$ conda update -c bioconda bedtools
``` 

g) deactivate the environment
``` {bash}  
$ conda deactivate
```
    
## 6 SLURM    
 
1. Create a new work folder and enter it    
    
``` {bash}
mkdir slurm
cd slurm/

```
    
2. Run squeue and squeue -u <username>. Are there any jobs submitted by you? NO 
    
``` {bash}
Run squeue and squeue -u Shilpa
```   
    
3. Submit jobs:
    
``` {bash}
for i in $(seq 1 15)
do
sbatch -J test_${i} -o output_${i}.out --wrap="sleep 20 && hostname"
done
```   
    
4. Check the job queue again
    
``` {bash}
squeue -u Shilpa
```   
    
    
5. When the job is done, compare the output to that of a local execution of hostname
    
``` {bash}
ls
cat output_14.out
$ hostname
```   
(Tells us that the tasks were executed on different nodes)
The node in the SLURM is bibigrid-worker-1-1-etqcgckjrdsvsxe
The node in local execution is bibigrid-master-etqcgckjrdsvsxe
    


    
    
    
    
    
    
## 7 Basics in R 
#### Data frames    

1. activate the core conda environment and start R session
```{r}
$ source /vol/COMPEPIWS/conda/miniconda3/bin/activate /vol/COMPEPIWS/conda/miniconda3/envs/core

$ R
```
    
2. Load the following libraries: “ggplot2”, “reshape2” and “GenomicRanges”
```{r} 
library(ggplot2)
library(reshape2)
library(GenomicRanges)
```
    
3. Read into R the following two tables: /vol/COMPEPIWS/groups/shared/*txt. Consider the first row as a header.
    
#go to the shared folder:
```{r}    
setwd("/vol/COMPEPIWS/groups/shared/")
```
#save the .txt files in a list:
```{r}    
filelist = list.files(pattern = ".txt")
```
#read:
```{r}    
table1 <- read.table("/vol/COMPEPIWS/groups/shared/kidney_14.5_mouse_1_2kbW_bed_counts.txt", header = TRUE)

table2 <- read.table("/vol/COMPEPIWS/groups/shared/liver_14.5_mouse_1_2kbW_bed_counts.txt", header = TRUE)
```
    
4. The dimension of each table :
    
```{r}
#dim of table1:
dim(table1) : 1362728 lines * 6 columns 
#dim of table2:
dim(table2) : 1362728 lines * 6 columns
```
    
5. Column names of each table :
    
```{r}
#col names table1:
colnames(table1) : "Chr"      "Start"    "End"      "H3K27me3" "H3K36me3" "H3K9me3"

#col names table2:
colnames(table2) : "Chr"      "Start"    "End"      "H3K27me3" "H3K36me3" "H3K9me3"
```
6. Genome length from each data set : 2000
    
```{r}
gen_len1 <- table1$End - table1$Start
gen_len2 <- table2$End - table2$Start
```

    
7. Create a data frame by concatenating vertically the two data sets. Add a new column called (cell_type) to annotate the rows of liver table with "liver" and the kidney’s rows with "kidney".
    
```{r}
#add new col in table1:
 table1["cell_type"] <- "kidney"
                        
```
                        
```{r}
                        
#add new col in table2:
 table2["cell_type"] <- "liver"
```
```{r}
#concatenate table1 and table2 vertically:
table <- rbind(table1, table2)
```
         
8. What is the dimension of the new data frame?
    
```{r}
dim(table): 2725456 * 7
```
         
9. Reshape the data frame (excluding chr, start, end) from wide into long format using "cell_type" as variable id. What are the dimensions of this dataframe?: 
Dimensions: 8176368 *  3

    
```{r}
table3 <- melt(table, id.vars = "cell_type", measure.vars = c("H3K27me3", "H3K36me3", "H3K9me3"))
```
```{r}
dim(table3)
```

10) Plot the density of each variable considering the cell_type in two different panels. Give different transparent color for each variable (hint: ggplot, geom_density, facet_wrap). Save the plot as a pdf (hint: ggsave).

```{r}
density_plot <- ggplot(table3, aes(x=enrichment, fill=mark)) + 
         geom_density(alpha=0.3) +
         facet_grid(cell_type~mark) +
         labs(x="Enrichment", y="Density")

ggsave(filename="density_plot.png", plot=p_density, device="png")
```
![](https://i.imgur.com/JbK0UxW.png)


11) The plots look very skewed, try to set the limit of X axis to 100 and re-plot. Save the plot as a pdf.

```{r}
scaled_density_plot <- ggplot(table3, aes(x=enrichment, fill=mark)) + 
                geom_density(alpha=0.3) +
                facet_grid(cell_type~mark) +
                labs(x="Enrichment", y="Density") +
                xlim(0, 100)

ggsave(filename="scaled_density_plot.png", plot=p_scaled_density, device="pdf")
```
![](https://i.imgur.com/PJPiVbh.png)

12) Create a new data frame by concatenating horizontally the two data sets. Consider only the following columns: H3K27me3, H3K36me3 and H3K9me3.

```{r}
table <- cbind(select(table1, c("H3K27me3", "H3K36me3", "H3K9me3")),
                   select(table2,  c("H3K27me3", "H3K36me3", "H3K9me3")) )
```
13) Rename the header of the new data frame into: H3K27me3_liver, H3K36me3_liver, ..........

```{r}
colnames(table) <- c("kidney_H3K27me3", "kidney_H3K36me3", "kidney_H3K9me3", "liver_H3K27me3", "liver_H3K36me3", "liver_H3K9me3")
```
14) Perform Principal Components Analysis on this data frame (hint: prcomp)
```{r}
PCA <- prcomp(table, scale=TRUE)
```
      
#### Genomic Ranges
1) Create two GenomicRange objects from the liver and kidney data frames using columns 4,5 and 6 as metadata
```{r}
gr_kidney <- GRanges(seqnames = table1$Chr,
                    ranges = IRanges(start=table1$Start, end=table1$End),
                    mcols= select(table1,c("H3K27me3", "H3K36me3", "H3K9me3")))

gr_liver <- GRanges(seqnames = table2$Chr,
                    ranges = IRanges(start=table2$Start, end=table2$End),
                    mcols= select(table2,c("H3K27me3", "H3K36me3", "H3K9me3")))
```

2) Calculate the total number of bases covered in these two objects
```{r}
coverage(gr_kidney)
coverage(gr_liver)
```
3) Subset one of the objects to have only chr2 lines
```{r}
gr_kidney[seqnames(gr_kidney) == "chr2"]
gr_liver[seqnames(gr_liver) == "chr2"]
```

4) Shift one of the object’s ranges by 100 upstream
```{r}
gr_kidney_shfit <- shift(gr_kidney, +100)
```      
    
5) Find overlapping regions between the two objects
```{r}
findOverlaps(gr_kidney_shfit, gr_liver)
```           
    
6) Make a GRangeList object from the two GRange objects
```{r}
GRangesList("Kidney" = gr_kidney, "Liver" = gr_liver)
``` 
            


# ChIP-seq Assay
## Nextflow ChIP-seq pipeline
### 1. Organize your raw data
a) Enter your working directory
``` {bash}
/vol/COMPEPIWS/groups/<group>/tasks
```   
b) List the raw fastq files under 
``` {bash}
/vol/COMPEPIWS/data/reduced/ChIP-seq
``` 

c) How many files per cell_type per time point per replicate are there? 
32 files per cell type 16 files per time point  2 files per replicate 
    
    
d) What are the histone marks you are dealing with?H3K27ac,H3K27me3,H3K4me1,H3K4me3,H3K9ac,H3K9me3,H3K36me3
    
e)  How many "control" per replicate do you have? : 2 or 4
    
f) Create a “data” folder in your group working directory and create symbolic links of the relevant files according to your group to the “data” folder as follow:
``` {bash}
ln -s  /vol/COMPEPIWS/data/reduced/ChIP-seq/*filename*  *symlinkname* 
```   
    
### 2. Generate the samplesheet
a) Enter the task folder
``` {bash}   
/vol/COMPEPIWS/groups/chipseq2/tasks
``` 
b) Prepare the samplesheet
```{bash}   
for i in /vol/COMPEPIWS/data/reduced/ChIP-seq/*liver*gz
do
group=$(basename ${i}|cut -f 1-3 -d _)
replicate=$(basename ${i}|cut -f 4 -d _)
fastq1=$i
fastq2=''
if [ "$(basename ${i}|cut -f 3 -d _)" == "control" ];
then
mark=''
control=''
else
mark=$(basename ${i}|cut -f 3 -d _)
control=$(basename ${i}|cut -f 1-2 -d _)"_control" 
fi
echo "$group,$replicate,$fastq1,,$mark,$control"
done > samplesheet.csv
```    
### 3. Running the pipeline
a) Based on the config file:
```
/vol/COMPEPIWS/pipelines/configs/chipseq.config.
```
i) Single end
ii) differential analysis,preseq,plot_profile,plot_fingerprint,spp
iii) Blacklist regions account for a significant portion of ChIP-seq reads, are driven by artifacts in genome assemblies, and removal of these regions is essential to removing noise in genomics assays.

Regions present : chr 18 & 19

iv) Tool used for mapping : BWA
    
b) Load Conda environment:
``` {bash}
source/vol/COMPEPIWS/conda/miniconda3/bin/activate/vol/COMPEPIWS/conda/miniconda3/envs/core
```
c) Run the pipeline
```{bash}
nextflow run nf-core/chipseq -r 1.2.2 -profile singularity -c /vol/COMPEPIWS/pipelines/configs/chipseq.config --input samplesheet.csv
```
4) To watch our jobs 
```
squeue --user=$ Shilpa
```
5) The output folder with the name "results"
```
/vol/COMPEPIWS/groups/chipseq2/results
```
6) Main steps of the pipeline for processing ChIP-seq : 
    Raw data,fastqc,Trimming,mapping (bwa),peak calling(MACS2)

## 2 Quality Control
2.1  Exploring the result folder
1.
```
/vol/COMPEPIWS/groups/chipseq2/results/pipeline_info/execution_report.html
```
a) Time taken to pipeline to finish : 2h 5m 35s
b) Step which took the most CPU/Memory usage : BWA_MEM (Mapping)
2. Aligned reads are located in :
```
/vol/COMPEPIWS/groups/chipseq2/results/bwa/mergedLibrary
```
3. Peak files are loacted in :
```
/vol/COMPEPIWS/groups/chipseq2/results/bwa/mergedLibrary/macs/narrowPeak
```
Tool used for Peak Calling: MACS2
4. Consensus peaks are located in :
```
/vol/COMPEPIWS/groups/chipseq2/results/bwa/mergedLibrary/macs/narrowPeak/consensus
```
5. Signal tracks are located in: 
```
/vol/COMPEPIWS/groups/chipseq2/results/bwa/mergedLibrary/bigwig
```

2.2 MultiQC
1. To explore diferent section of the MultiQC HTML report generated by the pipeline
```
vol/COMPEPIWS/groups/chipseq2/results/multiqc/narrowPeak
```
![](https://i.imgur.com/q4mWQKn.png)

    
2.
The quality score represents the probability that the corresponding nucleotide or base in the read is called incorrectly represented. It is calculated as Q = -10 x log10 P, where P is the probability that a base call is erroneous.From the figure below, note that there is a base 12 that has score around 23 which is lying below the green region. This indicates that base 12 of sample liver_14.5_H3K4me1_R2_T1 has poor quality, while rest of the bases of most of the samples are lying above green region and hence are of good quality.
    
![](https://i.imgur.com/63jaZdH.png)

We can see from the graph below that there are more number of reads for most of the samples that has average quality score more than 30. However, there are also fraction of reads of samples which are having average quality score is less than 30 (or lying in the red zone).

![](https://i.imgur.com/REbLT5w.png)

The figure below represents the graph of "number of reads" vs "sequence length". For instance we can see that liver_15.5_control_R2_T1 sample has more than 3M reads with sequence length of 50bp.
![](https://i.imgur.com/5MpYUM7.png)
    
Adapter contamination:
![](https://i.imgur.com/zFRgthR.png)

    
From given two figures below, we can note that the variable "number of reads" is changed. Note that, for the mentioned sample in figures, number of unique reads is increased while number of duplicate reads is decreased in "after trimming" as compared to "before trimming". This might be due to the fact that some of the duplicate reads after trimming become unique and hence number of unique reads is increased, though increase is not significant for the mentioned sample. Furthermore, we have lost some reads after trimming of around 100 of reads for the given sample in figures.
number of reads "before trimming (raw)":
![](https://i.imgur.com/6WuLcaM.png)

number of reads "after trimming":
![](https://i.imgur.com/WdXFte7.png)



3. Mapping efficiency: Mapping efficiency decreases with the increase of # unmapped reads as it give rise to low confidence on the score we calculate for further analysis. Note from the figure below, the sample liver_14.5_H3K4me3_R2 has high percentage of unmapped reads as compare to other samples (in terms of unmapped reads).Therefore this sample has comparatively low mapping efficiency.
![](https://i.imgur.com/7P6O3dI.png)


Duplication rate : If we have more duplicate reads then that will also affect our mapping quality. Large number of duplicates => high duplication rate. For instance, in the graph below, sample liver_14.5_H3K27ac_R2 has highest number of duplicate reads and therefore high duplication rate.
Q) what does these staright lines indicate in "Paired in sequencing" etc.?
![](https://i.imgur.com/2e3Ubac.png)


4. Do you get the similar number of peaks for all histone marks? : No. From the figure below, we can clearly compare that there are some marks with significantly larger number of peaks than other histone marks.
![](https://i.imgur.com/xZtNzdA.png)

    
   Overall: 
   Marks with the highest number of peaks : H3K4me1
   Marks with the lowest number of peaks : H3K27me3

5.What does FRiP score mean?
  It is generated by calculating the fraction of all mapped   reads that fall into the MACS2 called peak regions. A  read must overlap a peak by at least 20% to be counted. 
  
  Which marks get the highest/lowest FRiP scores? 
  ![](https://i.imgur.com/LlJlBCZ.png)
  Marks with highest FRiP scores : H3K4me3 (0.43)
  Marks with lowest FRiP scores : H3K4me1 (0.02)
  
  How is it related to the mark signal distribution (narrow/broad)?
  higher the frip score => broader the signal
    
    
2.3 deepTools2
1. Activate the core environment
```
source /vol/COMPEPIWS/conda/miniconda3/bin/activate /vol/COMPEPIWS/conda/miniconda3/envs/core
```
2. To Explore the batch effect between the samples
```
multiBigwigSummary bins -b liver_14.5_H3K27ac_R1.bigWig liver_14.5_H3K27ac_R2.bigWig liver_14.5_H3K27me3_R1.bigWig liver_14.5_H3K27me3_R2.bigWig liver_14.5_H3K36me3_R1.bigWig liver_14.5_H3K36me3_R2.bigWig liver_14.5_H3K4me1_R1.bigWig liver_14.5_H3K4me1_R2.bigWig liver_14.5_H3K4me3_R1.bigWig liver_14.5_H3K4me3_R2.bigWig liver_14.5_H3K9ac_R1.bigWig liver_14.5_H3K9ac_R2.bigWig liver_14.5_H3K9me3_R1.bigWig liver_14.5_H3K9me3_R2.bigWig liver_14.5_control_R1.bigWig liver_14.5_control_R2.bigWig liver_15.5_H3K27ac_R1.bigWig liver_15.5_H3K27ac_R2.bigWig liver_15.5_H3K27me3_R1.bigWig liver_15.5_H3K27me3_R2.bigWig liver_15.5_H3K36me3_R1.bigWig liver_15.5_H3K36me3_R2.bigWig liver_15.5_H3K4me1_R1.bigWig liver_15.5_H3K4me1_R2.bigWig liver_15.5_H3K4me3_R1.bigWig liver_15.5_H3K4me3_R2.bigWig liver_15.5_H3K9ac_R1.bigWig liver_15.5_H3K9ac_R2.bigWig liver_15.5_H3K9me3_R1.bigWig liver_15.5_H3K9me3_R2.bigWig liver_15.5_control_R1.bigWig liver_15.5_control_R2.bigWig -o chipseq2results.npz
```
number of bins found :15215
```{bash}
plotCorrelation -h
```
```{bash}
plotCorrelation -in chipseq2results.npz --whatToPlot heatmap --corMethod pearson -o heatmap.png
```
![](https://i.imgur.com/5SjRp3D.png)

```{bash}
plotPCA -in chipseq2results.npz -o pca.png
```
![](https://i.imgur.com/aa3ZnLV.png)

    
3. Marks that are similar to each other based on plotCorrelation outputs
H3K9ac ~ H3K4me3
H3K9ac ~ H3K4me1
H3K36me3 !~ H3K27me3

4. To separate signals (marks) from background: 
PATH > /vol/COMPEPIWS/groups/chipseq2/results/bwa/mergedLibrary/ to fingerprint.png
```{bash}
fingerprint: plotFingerprint -b liver_14.5_H3K9me3_R1.mLb.clN.sorted.bam liver_14.5_H3K9me3_R2.mLb.clN.sorted.bam liver_15.5_H3K9me3_R1.mLb.clN.sorted.bam liver_15.5_H3K9me3_R2.mLb.clN.sorted.bam -plot fingerprint_H3K9me3.png

plotFingerprint -b liver_14.5_H3K9me3_R1.mLb.clN.sorted.bam liver_14.5_H3K9me3_R2.mLb.clN.sorted.bam liver_15.5_H3K9me3_R1.mLb.clN.sorted.bam liver_15.5_H3K9me3_R2.mLb.clN.sorted.bam -plot fingerprint_H3K9me3.png
```
![](https://i.imgur.com/htBLBct.png)
    
![](https://i.imgur.com/QHNg3n9.png)
Not well separated from background for histone mark H3K9me3.

    
    
5. Marks enriched at TSS : H3K4me3, H3K9ac
    
```{bash}
computeMatrix scale-regions -S liver_14.5_H3K27ac_R1.bigWig liver_14.5_H3K27me3_R1.bigWig liver_14.5_H3K36me3_R1.bigWig liver_14.5_H3K4me1_R1.bigWig liver_14.5_H3K4me3_R1.bigWig liver_14.5_H3K9ac_R1.bigWig liver_14.5_H3K9me3_R1.bigWig -R /vol/COMPEPIWS/groups/chipseq2/results/genome/genome_genes.bed -a 1500 -b 1500 -o chipseq2matrix

plotProfile -m chipseq2matrix -o profileplot.png
```
    
![](https://i.imgur.com/xJAHCc7.png)

    
PATH > /vol/COMPEPIWS/groups/chipseq2/results/bwa/mergedLibrary/bigwig/ to profileplot.png and heatmap.png
    
## 3  Exploratory Analysis using IGV
1. IGV installation
```
/vol/COMPEPIWS/softwares/IGV_Linux_2.10.2/
```
2. To start IGV
```
/vol/COMPEPIWS/softwares/IGV_Linux_2.10.2/igv.sh
```
3. To open IGV session from : 
```
results/igv/narrowPeak/igv_session.xml
```
To remove consensus tracks : Right-click a track name and select Remove Tracks in the pop-up menu
    
4. To upload gene tracks
```
results/genome/genome_genes.bed
```
![](https://i.imgur.com/XUMRpND.png)
    
5. Do you see clear peaks for all marks? : Yes

   Sharper marks : H3K9ac,H3K4me3,H3K4me1,H3K36me3,H3K27me3,H3K27ac

    ![](https://i.imgur.com/yfYH3q9.png)

    
    
    
   Broader marks : H3K9me3 & Control

    ![](https://i.imgur.com/y4qz4yI.png)
6. Coloring the tracks as we want (e.g. per mark). Right click on the track name.
   
    
   ![](https://i.imgur.com/QvETQSV.png)
7. Marks that overlap with TSS : H3K27ac, H3K4me3, H3K9ac
   Marks that overlap with gene body : H3K9ac, H3K4me1, H3K36me3, H3K27me3
  
      ![](https://i.imgur.com/CakzSAJ.png)



8. Do you always see an agreement between the called peaks of the two replicates of the same mark in the same cell type? What do you understand from this plot (choose two of the sharp marks)

No, this is always not the case. As we can see in the figure below that for the mark H3K27ac R1 and R2 are not agreeing on the same region. This may arise due to the signal noise or the may be due to contamination of the solution while preparing the sample.
    
![](https://i.imgur.com/LUfoX9K.png)
![](https://i.imgur.com/xwqTpAf.png)


    
    
## 4 Chromatin segmentation with ChromHMM

4.1 ChromHMM Installation
```
/vol/COMPEPIWS/softwares/ChromHMM/
```
    
4.2 Running ChromHMM
1) Activate Core conda environment
```
source /vol/COMPEPIWS/conda/miniconda3/bin/activate /vol/COMPEPIWS/conda/miniconda3/envs/core
```
2) Create a folder called segmentation under your work group folder, and enter it
```
mkdir segmentation
```

3) Create a folder called “inputs” 
```
mkdir inputs
```
Create symbolic links of all aligned read files to this folder (there should be 64 bam files).
```
for i in /vol/COMPEPIWS/groups/chipseq2/results/bwa/mergedLibrary/*liver*bam; do ln -s $i /vol/COMPEPIWS/groups/chipseq2/segmentation/inputs/ ; done
```
4) Prepare a tab delimited file “cellmarkfiletable.txt”
```
for i in *liver*bam *kidney*bam 
do
cell_type=$(basename ${i}|cut -f 1-2 -d _)
histone_mark=$(basename ${i}|cut -f 3 -d _)
celltype_stage=$(basename ${i}|cut -f 1-4 -d _)
if [ "$(basename ${i}|cut -f 1-2 -d _)" == "liver_14.5" ];
then
control='liver_14.5_control_R1.mLb.clN.sorted.bam'
elif [ "$(basename ${i}|cut -f 1-2 -d _)" == "liver_15.5" ];
then
control='liver_15.5_control_R1.mLb.clN.sorted.bam'
elif [ "$(basename ${i}|cut -f 1-2 -d _)" == "kidney_14.5" ];
then
control='kidney_14.5_control_R1.mLb.clN.sorted.bam'
elif [ "$(basename ${i}|cut -f 1-2 -d _)" == "kidney_15.5" ];
then
control='kidney_15.5_control_R1.mLb.clN.sorted.bam'
fi
printf '%s\t%s\t%s\t%s\n' "$cell_type" "$histone_mark" "$celltype_stage" "$control"
done > cellmarkfiletable.txt
```

5)  Run the BinarizeBam command with 200bp as a bin size
```
java -Xmx3000M -Djava.awt.headless=true -jar /vol/COMPEPIWS/softwares/ChromHMM/ChromHMM.jar BinarizeBam -b 200 /vol/COMPEPIWS/softwares/ChromHMM/CHROMSIZES/mm10.txt /vol/COMPEPIWS/groups/chipseq2/segmentation/inputs inputs/cellmarkfiletable.txt outputs
```
6) Run the LearnModel command with different number of states: 5-15
```
for i in $(seq 5 15)
do
sbatch --job-name=learnmodel -o output_${i}.out --wrap="java -Xmx3000M -Djava.awt.headless=true -jar /vol/COMPEPIWS/softwares/ChromHMM/ChromHMM.jar LearnModel /vol/COMPEPIWS/groups/chipseq2/segmentation/outputs /vol/COMPEPIWS/groups/chipseq2/segmentation/learnmodel_out $i mm10"
done
```
7) Explore the output files called webpage_$state.html. Especially, emission heatmap, Ref Seq_TSS_neighborhood and overlap figures. Comment on the “new” states you get.
    
![](https://i.imgur.com/CEFKOOH.png)

![](https://i.imgur.com/cfCQELY.png)

![](https://i.imgur.com/4sH7lyN.png)


    
8) Compare the different models to assess what is the proper minimum number of states: 10 number of states are appropriate for
```
java -Xmx3000M -Djava.awt.headless=true -jar /vol/COMPEPIWS/softwares/ChromHMM/ChromHMM.jar CompareModels emissions_15.txt /vol/COMPEPIWS/groups/chipseq2/segmentation/learnmodel_out/ compare_model
```
9) Repeat step 5 to merge the stages of each cell type in order to produce only 2 segmentation files (1 for each cell type).
```
java -Xmx3000M -Djava.awt.headless=true -jar /vol/COMPEPIWS/softwares/ChromHMM/ChromHMM.jar BinarizeBam -b 200 /vol/COMPEPIWS/softwares/ChromHMM/CHROMSIZES/mm10.txt /vol/COMPEPIWS/groups/chipseq2/segmentation/inputs segmentation/cellmarkfiletable.txt merge_outputs
```
10) Repeat steps 6 with the new “cellmarkfiletable” for only 15 states model
```
sbatch --job-name=learnmodel -o output_${i}.out --wrap="java -Xmx3000M -Djava.awt.headless=true -jar /vol/COMPEPIWS/softwares/ChromHMM/ChromHMM.jar LearnModel /vol/COMPEPIWS/groups/chipseq2/merge_outputs /vol/COMPEPIWS/groups/chipseq2/segmentation/learnmodel_merge_out 15 mm10"
```
11) Based on the emission heatmap, overlap and enrichment of this new 15-state model, assigning biological
labels to the states 
![](https://i.imgur.com/KMdp99B.png)
![](https://i.imgur.com/9ZaY8Fg.png)




12) Rename the states in column 4 of the files *dense.bed (or *segments.bed) using the biological states from the previous question (consider the command MakeBrowserFiles from ChromHMM or awk in bash to do the task). Note: Redirect the output into new files in a subfolder called “relabeled” to avoid overwriting the original files.
For liver:
```
java -Xmx3000M -Djava.awt.headless=true -jar /vol/COMPEPIWS/softwares/ChromHMM/ChromHMM.jar MakeBrowserFiles -m /vol/COMPEPIWS/groups/chipseq2/segmentation/seg_name.txt /vol/COMPEPIWS/groups/chipseq2/segmentation/learnmodel_merge_out/liver_15_segments.bed seg_chipseq2 /vol/COMPEPIWS/groups/chipseq2/segmentation/relabelled/liver_relabelled
```
For kidney : 
```
java -Xmx3000M -Djava.awt.headless=true -jar /vol/COMPEPIWS/softwares/ChromHMM/ChromHMM.jar MakeBrowserFiles -m /vol/COMPEPIWS/groups/chipseq2/segmentation/seg_name.txt /vol/COMPEPIWS/groups/chipseq2/segmentation/learnmodel_merge_out/kidney_15_segments.bed seg_chipseq2 /vol/COMPEPIWS/groups/chipseq2/segmentation/relabelled/kidney_relabelled
```
13) Generating the emission heatmap with the new labels using the command "reorder"
```
java -Xmx3000M -Djava.awt.headless=true -jar /vol/COMPEPIWS/softwares/ChromHMM/ChromHMM.jar Reorder -m /vol/COMPEPIWS/groups/chipseq2/segmentation/seg_name.txt /vol/COMPEPIWS/groups/chipseq2/segmentation/learnmodel_merge_out/model_15.txt /vol/COMPEPIWS/groups/chipseq2/segmentation/relabelled/output  
```
14) Copy the relabeled segmentation files of Q12 and the emission heatmap of Q13 into the shared folder under
your group sub-directory /vol/COMPEPIWS/groups/shared/<group>/segmentation
```
for i in /vol/COMPEPIWS/groups/chipseq2/segmentation/relabelled/*.bed; do ln -s $i /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/segmentation; done
    
for i in /vol/COMPEPIWS/groups/chipseq2/segmentation/relabelled/output/*emissions*png; do ln -s $i /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/segmentation; done
```


## 5 Exploration of chromatin state segmentation using IGV
1. Which state(s) is(are) more likely overlapping with TSS?
![](https://i.imgur.com/2cUWWpZ.png)
    
2. Some differential staes in kidney and liver samples.
![](https://i.imgur.com/4tni12B.png)
    
## 6 Exploratory Analysis using R
1) Activate the “core” environment.
```
source /vol/COMPEPIWS/conda/miniconda3/bin/activate /vol/COMPEPIWS/conda/miniconda3/envs/core
Start R
```
2) Start R
3) Read all *(reordered)dense.bed files from task 4.2.12 from “segmentation” directory results into R. (hint:
use list.files, lapply, read.table functions)
```
dir <-"vol/COMPEPIWS/groups/chipseq2/segmentation/relabelled/"
regex <- "*dense.bed$"
filenames <- list.files(path=dir,
                        pattern=regex,
                        full.names=TRUE)
filenames
BEDs <- lapply(filenames,
               read.table,
               header=FALSE,
               sep="\t",
               stringsAsFactors=FALSE,
               skip=1)

BEDs[[1]]$sample <- "kidney"
BEDs[[2]]$sample <- "liver"
```
4) Merge the results into one data frame (hint: rbind function)
```
df <- rbind(BEDs[[1]], BEDs[[2]])
```

# Integrative analysis
## 2. Integrative data exploration in IGV 
1) What are the methylation states of the promoters and the gene body in general (high, low)?
   In the figure below (zoomed in), we can see that the region of interest (genomic region: chr18:31,755,635-31,775,374) has an unmethylated promoter  for both kidney and liver sequences. While the gene body corresponding to this promoter region is highly methylated (HMR)  to maintain the transcription process and hinder the binding of other proteins to the gene body. Moreover, gene body represents strong transcription (Tx_S).

![](https://i.imgur.com/tDAjegn.png)
    
For an other instance (figure given below), it is clearly visible for the full length of gene 'Slc25a46' that promoters regions for kidney and liver samples are mostly unmethylated while their gene body lie under HMR.
![](https://i.imgur.com/9Tohbh6.png)
    
So in general, we can conclude that the promoter regions lie under UMR while their corresponding gene body lies under HMR to maintain the transcription process and complete it successfully. 
Note: Regions representing Tx_S as well as Tx_W are both HMR so we observe that level of trancription does not affect the methylation state of gene body.




2) Which states from DNA methylation segmentation (MethylSeekR) overlap with chromatin segmentation (ChromHMM)?
In the figure below, we can clearly see that for liver segmentation, state NS (no signal) is overlapping with HMR, Pr_A (active promoter) region is overlapping with UMR, Tx_W as well as Tx_S are overlapping with HMR. Kidney segmentation also follows the similar overlapping patterns.
![](https://i.imgur.com/dVITP5d.png)
    
Note from the figure below, we observe that there is a state 14_Het_H3K9me3 (H3K9me3 heterochromatin) that overlaps with HMR in liver sample while in kidney sample 14_Het_H3K9me3 overlaps with PMD. 
![](https://i.imgur.com/dtrkY8v.png)
    
Furthermore, there is a small region in the kidney sample with state NS that also overlaps with LMR. 
![](https://i.imgur.com/0p4qmUo.png)


3) Do you observe overlap between the accessibility peaks and specific chromatin states? Which ones?
Yes, there are  accessibility peaks that overlaps with specific chromatin states. As we can see from the figure below, accessibility peaks are overlapping with Pr_A (active promoter) state. 
![](https://i.imgur.com/MsSIHiI.png)
    
Furthermore, we have a region of interest in the liver genome where accessibility peak is overlapping with Pr_B (bivalent promoter: Histone H3 lysine 4 monomethylation [H3K4me1] or trimethylation [H3K4me3] and repressing (e.g. H3K27me3)). So here we have different states of chromatin in liver and kidney genome that overlap with accessibility peaks.
![](https://i.imgur.com/3683eJE.png)

Besides from accessibility peak - overlapping - Pr_A, we encounter the genomic region in liver sample where the chromatin state is SEnh_pxS (Strong,TSS -prox Enhancer)	that overlaps with accessibility peak.
For comparison purpose, in kidney sample, peaks also appears that overlaps with Tx_W (weak transcription) but not as sharp as liver sample.
![](https://i.imgur.com/gnDWwfx.png)
    
So we can say that mostly active promoter regions are overlaps with accessibility peak because they are associated with chromtin accessibility. While there may also (may be few) exist other states than promoter regions that overlaps with accessibility peak. 

CONSENSUS PEAKS: Peaks where both kindey and liver sample have chromatin accessibility.




4) Can you find some examples for unexpressed genes. What are the methylation states of their promoters and
bodies? Which histone marks are enriched at the promoters and bodies of these genes?
    
Yes, there are some gene which are not expressed. 
(a) Gene Gm3734: From the figure below, we see that the gene Gm3734 is not expressed as there are no signals in from RNAseq for both liver and kidney. Methylation states of its promoter and gene is PMD (partially methylated domain). 
![](https://i.imgur.com/SroK1Tj.png)


(b) Gene Hdhd1a: We found another gene which is unexpressed. This gene (promoter as well as gene body) is lying in PMD. Furthermore, we can see for both liver and kidney sample, that histone marks H3k27me3 and H3k9me3 show higher signals on the gene body and also on promoter regions.
![](https://i.imgur.com/fJL5upL.png)

    
c) Gene Ftmt: Figure below explores an another instance of unexpressed gene with no RNAseq signals. Further, the promoter and gene body of the region are lying under HMR (for both kidney and liver samples) (from WGBS signals and segmentation). The histone marks that we observe in this region are H3K27me3 and H3K9me3 (from chipseq signals). 
![](https://i.imgur.com/caMHTOK.png)

5) DNA & chromatin accessibility = negatively associated
DNA & gene expression = negatively associated

![](https://i.imgur.com/ph8hkfm.png)
    
    
![](https://i.imgur.com/t25CjFY.png)

6) Do you observe positive or negative associations between DNA methylation, chromatin accessibility, individual histone marks and gene expression? Do the association differ in certain types of regions (e.g. promoters, gene bodies, enhancers)?
Negative association between DNA methylation and the marked histone marks (figure below).
![](https://i.imgur.com/xMwob8b.png)

Negative association between DNA methylation and gene expression (figure below).
![](https://i.imgur.com/p5KVLQ5.png)
    
Positive association between DNA methylation (H3K4me3) and gene expression.
![](https://i.imgur.com/Un0kfgr.png)





## 3. Integrative analysis
1) Summarize signals over regions of interest
a) Define regions of interest (ROI) – generate BED3 files1 for one or more of the following regions:
    ROI : promoter region (Pr)
```
grep "Pr" /vol/COMPEPIWS/groups/chipseq2/segmentation/relabelled/liver_relabelled_dense.bed | awk '{print $0}' > /vol/COMPEPIWS/groups/chipseq2/analysis/Pmr.bed
    
computeMatrix scale-regions -S /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/signals/liver_14.5_H3K4me3_R1.bigWig /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq3/signals/kidney_14.5_H3K9me3_R1.bigWig -R Pmr.bed -o Pr_Comp_Enh
    
```    
b) Combine computeMatrix with plotHeatmap/plotProfile from deepTools22 to generate plots. e.g. histone mark/accessibility/DNA methylation signals across all/some of the above mentioned ROI. What do you observe? Explain your findings from a gene-regulation perspective.
    

```
plotProfile --matrixFile Pr_Comp_Enh --outFileName plotprofile_Pmr_CompEnh.png

plotHeatmap --matrixFile Pr_Comp_Enh --outFileName plotheatmap_Comp_Pmr_Enh.png
    
```   
    
The figure below predicts that the histone mark H3K4me3 gives rise to gene expression while the histone mark H3K9me3 is not associated with gene expression. The black spots denotes the noise that should be rectified. Blue region denotes the gene expression, red region denote repressed region.
![](https://i.imgur.com/dYw6ppi.png)


    
2) Chromatin states and DNA methylation
a) Average DNA-methylation in the chromatin states (15-state model) of the corresponding samples
    
```
bedtools intersect -wa -wb -a "/vol/COMPEPIWS/groups/shared/WGBS/wgbs1/signal/rnbeads_sample_01_kidney.bedGraph" -b "/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq3/segmentation/kidney_relabeledstate15_dense.bed_dense.bed"    

bedtools groupby -i liver_avg_DNAme.bed -g 8 -c 4 -o mean > liver_group.bed

In R:
ggplot(obj, aes(x=Histone, y=liver_me)) + geom_boxplot(fill="pink")
```
b) 
![](https://i.imgur.com/QxdWtJi.png)
![](https://i.imgur.com/HC5xD2E.png)
    
c) From the figure we can see that the state Pr_A (active promoter) is low methylated in both kidney and liver with histone marks: H3K27ac, H3K9ac, H3K4me3. While transcription states are highly methylated with histone marks: H3K36me3.
    
d) Yes, there is a difference in the average methylation between two cell types. For instance, Pr_A is more methylated in kidney sample as compared to liver sample. Furthermore,Enh_dsS and Enh_pxS are more methylated in liver sample as compared to kidney sample.


4) Differentially Methylated Regions (DMRs) between kidney and liver: overlap with chromatin states and
accessibility
a)  

```  
cut -f 16-19 "/vol/COMPEPIWS/groups/shared/WGBS/wgbs1/differential/methdiff.tsv" |sed 's/"//g'  > mythl.bed

grep Het "/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq2/segmentation/kidney_relabelled_dense.bed" > kidney_DMR.bed

bedtools intersect -wa -wb -a mythl.bed -b kidney_DMR.bed | wc -l
 ``` 
528 and 1061 (for liver and kidney respectively).
    

![](https://i.imgur.com/3mhiJgh.png)

    
