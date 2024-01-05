# Document Outline     
Herein, I will seek to add some additional support to the potential lncRNAs
within our dataset. By intersecting the transcripts with GRO-seq data, we will
provide additional support for nascent transcription at the locus. By
overlaying with p300+H3K4me1 data we will work to filter out those that may be
enhancer related, and by looking for conservation blocks and elevated PhyloP
signals, we may prioritise those functional-putative genes amongst our set.     


# Analysis

## General Charactersitics of our set

Average read length, min and max transcript length. 
```
less -S lincRNA-non-coding-set-superfile2.tab | awk '{ sum += $4 } END { if (NR > 0) print sum / NR }' 
less -S lincRNA-non-coding-set-superfile2.tab | sort -k4,4 -rn | head/tail
```
* 1856bp avg. length, 212bp min, 10141 max (monoexonic too).      
* Median of 3904

Average coverage level.    
`less -S lincRNA-non-coding-set-superfile2.tab | awk '{ sum += $2 } END { if (NR > 0) print sum / NR }' | head`
* 5.3 transcripts - huge skew by the 795x transcript, the vast majority are below 20x.            



## GRO-seq intersect      


Intersect with GRO-seq data. More information regarding the GRO-seq datasets used can be found in Pseudogene-GROseq-additional-support.markdown        
* First, remove pseudogene IDs, and, create a lncRNA .bed file.    

```     
less -S $PSSET | awk '{print $11}' | sed 's/["";]//g' | fgrep -w -v -f - lincRNA-non-coding-set-superfile.tab > lincRNA-non-coding-set-superfile2.tab   

awk '{print $1}' lincRNA-non-coding-set-superfile2.tab | fgrep -w -f - ../../../true-unknown-noncoding.bed > lincRNA-non-coding-set.bed

awk '{print $1}' lincRNA-non-coding-set-superfile2.tab | fgrep -w -f -
../../../true-unknown-noncoding.bed | awk '$8 == "transcript" {print $0}' >
lincRNA-non-coding-set-transcripts.bed     
```     
*631 transcripts 

Now intersect with the GRO-seq data.   

```
sort -k1,1 -k2,2n lincRNA-non-coding-set-transcripts.bed | bedtools intersect
-a - -b $GRO | awk '{print $11}' | sed 's/["";]//g' | sort | uniq | wc -l

sort -k1,1 -k2,2n lincRNA-non-coding-set-transcripts.bed | bedtools intersect -a - -b $GRO2 | awk '{print $11}' | sed 's/["";]//g' | sort | uniq | wc -l
```
* 219 transcripts are supported by GRO-seq dataset 1. 238 by GRO-seq dataset 2. GSE27037.         


## Enhancer intersection     

Greenberg 2010 enhancer datasets - p300 and H3K4me1 Chip-Seq. Filter out
lncRNAs that may be enhancer based, and perhaps confirm those with NO-CRE
signature, but nonetheless resemble enhancer transcripts.  


Can filter out bedgraph files by selecting only those regions with binding scores above baseline - perhaps something over 6?    
`less -S GSM530203_H3K4Me1_H3K4Me1_un_B2_E120.bedgraph | awk '$4 > 6 {print $0}'`       

Intersect with CPB, and K3K4me1 in both untreated cells and KCL stimulated.    
`bedtools intersect -sorted -a $lncbed -b GSM530204_H3K4Me1_H3K4Me1_KCl_B2_E120.bed -wb` 
* Iterted over all the files as above.     
* No binding events over baseline levels ~1 and under. Peak reads at ~150+ 



## PhyloP conservation signature    

One valuable method which can be used to parse the 600 or so lincRNA candidates
is to look for the transcripts with elevated conservation (PhyloP) levels.
Evolution, by virtue, should keep those elements which are necessary.    

I will download the PhyloP 60way from UCSC and intersect it with the lincRNA candidates.    

Download from UCSC mm10 golden path (60 way)     

`/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/UCSC-PhyloP-Cons/PhyloP/`     

Convert from bigwig to bedgraph, and then into an intersectable bed format. 

```     
# Use Kent Tool utility bigwigtoBedGraph 

# From bedgraph to bed with this simple awk script    
awk '{print $1 "\t" $2 "\t" $3 "\t" $4}')    
```

Intersect with PhyloP .bed.     
```   
lncbed=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-tru
ly/non-coding-set/NO-CRE-support/Other-Genetic-Elements/lincRNA-ncRNA/lincRNA-non-coding-set.bed    

bedtools intersect -sorted -a $lncbed -b ./mm10.60way.phyloP60way.bed -wb > tmp.out
```     

Split transcripts with phyloP scores over three into their own subset. This can be viewed as the highest confidence transcripts. 
* A large question is whether within these, we have some structural RNAs?     i
* How to reconcile high phyloP but very low expression? 
`less lincRNA-set-NC-PhyloP.bed |  awk '$19 >= 3 {print $0}' > lincRNA-set-NC-PhyloP-Over3.bed`    
* 151 transcripts within this set.      

Create a 'SuperFile' for these too.   
`less -S lincRNA-set-NC-PhyloP-Over3.bed | awk '{print $11}' | sort | uniq |
sed 's/["";]//g' | fgrep -w -f - ../lincRNA-non-coding-set-superfile2.tab >
lincRNA-set-NC-PhyloP-Over3-Superfile.tab`      

Script to easily search transcript IDs with the highest sequence conservation.     
less -S tmp.out | awk '{print $11, $19}' | sort -k2,2 -rn | awk '$2 > 3 {print $0}' | cut -d" " -f1 | sort | uniq | less -S


It is clear that parsing these is a very enjoyable and potentially fruitful
endevour. Some elements have truly astonishing levels of conservation would
should be addressed.    


A few examples to visualise with UCSC
STRG.14535.1 - Conserved exon-intron structure, open atac peaks.     
STRG.197.1 - High conservation 
STRG.21637.2 - High cons, 20x coverage..   







