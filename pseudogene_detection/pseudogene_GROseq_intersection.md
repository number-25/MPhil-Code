# Document Outline 
Herein, I will intersect the novel pseudogenes in our dataset with GRO-seq
datasets in order to look for additional lines of transcriptional support for
these transcripts, and thus bolster claims that they are in fact pseudogenes.   

There are two main public GRO-seq datasets that I will pull from. 

The first is (GSE27037)[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27037] as
described in this publication ("Regulating RNA polymerase pausing and
transcription elongation in embryonic stem cells
")[https://pubmed.ncbi.nlm.nih.gov/21460038/]. 

The second is (GSE43390)[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43390] as
described in this publication ("Pausing of RNA polymerase II regulates
mammalian developmental potential through control of signaling
networks")[https://pubmed.ncbi.nlm.nih.gov/25773599/#:~:text=of%20signaling%20networks-,Pausing%20of%20RNA%20polymerase%20II%20regulates%20mammalian%20developmental%20potential%20through,Mol%20Cell.]          

## Analysis 

Protein-coding + Non-coding candidates intersection. 

First do GSE27037.       
```
GRO=/media/labpc/Disk-2/Bioinformatics-Computational/GEO-Datasets/GRO-seq/Regulating-RNA-pol-Pausing-GSE27037/GSE27037_MEF_GROseq-sorted.bed

GRO2=/media/labpc/Disk-2/Bioinformatics-Computational/GEO-Datasets/GRO-seq/Regulating-RNA-pol-Pausing-GSE27037/GSE27037_MESC_GROseq-sorted.bed

PSSET=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Pseudogene-Pred/NC-PC-Pseudogene-set-transcript.bed   

bedtools intersect -a $PSSET -b $GRO -wb > PS-GRO-seq-intersect-GSE27037_MEF.bed     

bedtools intersect -a $PSSET -b $GRO2 -wb > PS-GRO-seq-intersect-GSE27037_MESC.bed

bedtools intersect -a $PSSET -b $GRO -wb | awk '{print $11}' | sed 's/["";]//g' | sort | uniq | wc -l       
bedtools intersect -a $PSSET -b $GRO2 -wb | awk '{print $11}' | sed 's/["";]//g' | sort | uniq | wc -l        
```      
* 15 pseudogenes are supported by GRO-seq dataset 1, 29 transcripts supported by GRO-seq dataset 2.        


Now do GSE43390.      

```    
G1F=/media/labpc/Disk-2/Bioinformatics-Computational/GEO-Datasets/GRO-seq/Pausing-RNA-pol-II-GSE43390/GSE43390_Adelman_Bl6_mESC_C2cells_2i_4OHT_startRNA-seq_5pr_allReps_norm_forward.bed

G1R=/media/labpc/Disk-2/Bioinformatics-Computational/GEO-Datasets/GRO-seq/Pausing-RNA-pol-II-GSE43390/GSE43390_Adelman_Bl6_mESC_C2cells_2i_4OHT_startRNA-seq_5pr_allReps_norm_reverse.bed    

G2F=/media/labpc/Disk-2/Bioinformatics-Computational/GEO-Datasets/GRO-seq/Pausing-RNA-pol-II-GSE43390/GSE43390_Adelman_Bl6_mESC_C2cells_2i_Ctrl_startRNA-seq_5pr_allReps_norm_forward.bed    

G2R=/media/labpc/Disk-2/Bioinformatics-Computational/GEO-Datasets/GRO-seq/Pausing-RNA-pol-II-GSE43390/GSE43390_Adelman_Bl6_mESC_C2cells_2i_Ctrl_startRNA-seq_5pr_allReps_norm_reverse.bed     

bedtools intersect -a $PSSET -b $G1F -wb > PS-GRO-seq-intersect-GSE43390_G1F.bed    

bedtools intersect -a $PSSET -b $G1R -wb > PS-GRO-seq-intersect-GSE43390_G1R.bed     

bedtools intersect -a $PSSET -b $G2F -wb > PS-GRO-seq-intersect-GSE43390_G2F.bed     

bedtools intersect -a $PSSET -b $G2R -wb > PS-GRO-seq-intersect-GSE43390_G2R.bed     
```     


List of transcripts which have GRO-seq overlap across datasets. Ideally, for visualisation
purposes, I'd like to find a pseudogene which is covered by short-reads, and
more than 1 GRO-seq dataset.     


GSE27037-MEF  
STRG.3249.2 chr10   85818326  MESC also
STRG.3955.1 chr11   8338063   MESC also
STRG.4019.1 chr11   20114602  MESC, GSE43390_G1F, G2F, Short-read also (made graph)     
STRG.8189.1 chr13   55389334  MESC, short-read also
STRG.12346.1 chr16   46866027 MESC also
STRG.17463.1 chr2    90355301 Short read also
STRG.26035.1 chr6   92084828  MESC, GSE43390_G1F, G2F, Short-read also     (made graph)    

As above - MESC
STRG.12737.1  chr17   6492267 Short read also
STRG.17533.1  chr2    17941088  G2R, short read also,
STRG.30854  chr9    6817359 

GSE43390_G1F     
STRG.17533.1  chr2    92812434  G1R, G2F also
STRG.28365.1  chr7    105252340 G1R, G2F also,

GSE43390_G1R     
STRG.14624.1  chr18   34532830  G2R also, 

G2F
STRG.16503.1  chr2    17941006
STRG.18562.1  chr2    162919333 Short read also

3'UTR Scan Predicted Pseudogene intersection.     
* These will be perhaps the most valuable for these pseudogenes were not discovered with transcriptional evidence, but computationally instead.    

```    
UTRPS=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Pseudogene-Pred/3UTR-scan/high-conf-unique/PC+NC-predicted-PS-manual-complete-sorted.bed     
```
