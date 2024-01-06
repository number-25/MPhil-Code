# Document Outline

Herein, I will attempt to perform a second pseudogene verification strategy,
based upon the observations within the directRNA dataset (noncoding blastn
paired with manual BLAT). I will be addressing the question which asks whether
the pseudogenes detected in the dRNA dataset, have not been detected previously
due to inadaquencies of short read sequencing and its transcript
reconstruction. 

I pondered whether I could simulate short-read data from the pseudogenes, and
their parents genes annnotations, and map them to the genome, looking at
whether the short reads are able to distinguish the novel pseudogenes from
their parents, and the effects of read-depth on their reconstruction and thus
discovery. The read depth component is particularly interesting as many of
the novel pseudogenes have very, very low amounts of coverage in our dataset.   

* Is discovery of novel pseudogenes a function of read length, and read depth?    
* Longer reads carry more information and thus more discriminatory power - higher resolution.    


### Analysis

Only one RNAseq simulation software will be utilized, as all others are either dated, incomplete, or dedicated to genome sequencing data.         

1. Polyester by Jeff Leek's (group)[https://academic.oup.com/bioinformatics/article/31/17/2778/183245?login=true].     

Runners up:
RNASeqReadSimulator: (A Simple RNA-Seq Read Simulator)[http://alumni.cs.ucr.edu/~liw/rnaseqreadsimulator.html] 
- This is probably the most dated of the three, it runs on Python2.7 and
strictly requires BED12 files as input. Nonethless, it has some simply and
clever features.     
- An R program - elaborate and precise, with many options - likely to be the most well rounded and consistent option.     
Random reads from the (BBMap
suite)[https://github.com/BioInfoTools/BBMap/blob/master/sh/randomreads.sh].) -
Appear to be mainly used for genomic reads, and it doesn't preserve exon-intron
structure.      


## Analysis     

Create a custom .gtf file containing the transcript structure of the novel
pseudogenic transcripts, as well as their parent genes. This file will be
supplied to Polyester along with a .fasta file in order to generate the
simulated reads.     

Non-coding set parent pseudogenes.      
```     
PSLST=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/lo
ng-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/non-coding-set/pseudogenes/nc-pseudogenes-hiconf.lst    

vm24=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm24/gencode.vM24.annotation.gtf      

cut -f2 nc-pseudogenes-hiconf.lst | fgrep -w -f - $vm24 > custom-NC-set-PS-parents.gtf 

cut -f2 $PSLST | fgrep -w -f - $vm24 | awk '$3 == "gene" {print $0}' | awk '{print $14}' | sort | uniq | grep -v "-" | wc -l     
```      
* Confirmed that the 50 non-coding set pseudogene 'parents' are within the custom.gtf.          
* Within the custom annotation are isoform specific annotations also e.g. Gs2-201, Gs2-202 etc.      

Protein-coding set parents pseudogenes.     

*Using vm25 here for more completeness.    
```    
PCLST=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-An
notation/Trinotate/TransDecoder/blastn/pseudogenes/pc-candidates-pseudogenes-blastn.lst     

vm25=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm25/vm25-gencode.gtf     

cut -f2 pc-candidates-pseudogenes-blastn.lst | fgrep -w -f - $vm25 | awk '$3 ==
"gene" {print $0}' | awk '{print $14}' | sort | uniq | grep -v "-" | wc -l     

cut -f2 $PCLST | fgrep -w -f - $vm25 > custom-PC-set-PS-parents.gtf
```   
* Confirmed that the 12 protein-coding set pseudogene 'parents' are within the
custom.gtf.     
* Within these annotations are also known pseudogenes of the functional parent
gene, which in my opinion is of major benefit, as it mimics the natural
situation in which multiple, highly similar transcripts are 'competing' for
primary mapping.      


Merge both newly created cutom .gtf files.     
`cat custom-NC-set-PS-parents.gtf custom-PC-set-PS-parents.gtf > custom-PC-NC-PS-parents.gtf`     
* All 62 parents genes are here.     


The next step is to create a .gtf file for the unannotated pseudogenes within our dRNA dataset.    

Convert .bed files to .gtf files, then merge both files.     
```     
awk '{print $1"\t"$7"\t"$8"\t"($2+1)"\t"$3"\t"$5"\t"$6"\t"$9"\t"(substr($0, index($0,$10)))}' nc-pseudogenes-hiconf.bed > nc-pseudogenes-hiconf.gtf    

NCPSGTF=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/non-coding-set/pseudogenes/nc-pseudogenes-hiconf.gtf     

awk '{print $1"\t"$7"\t"$8"\t"($2+1)"\t"$3"\t"$5"\t"$6"\t"$9"\t"(substr($0,
index($0,$10)))}' pc-candidates-pseudogenes-blastn.bed >
pc-candidates-pseudogenes-blastn.gtf     

PCPSGTF=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/TransDecoder/blastn/pseudogenes/pc-candidates-pseudogenes-blastn.gtf     

cat $NCPSGTF $PCPSGTF > custom-PC-NC-novel-set.gtf    
```

Create final .gtf for input to polyester, as well as a corresponding .fasta file.     

`cat custom-PC-NC-PS-parents.gtf custom-PC-NC-novel-set.gtf > custom-polyester-input-ALL.gtf`       

* Make sure to include exon co-ordinates in fasta header, and ensure that the spliced transcript nucleotides are output.     
`~/Bioinformatics-Programs/gffread/gffread -W -w custom-polyester-input-ALL-splicedtranscripts.fa -g $mm10 custom-polyester-input-ALL.gtf`    


### R - polyester package.    

A few iterations of transcript simulation will be performed, with progressively
increasing levels of coverage, until all transcripts are able to be
reconstructed/mapped. By assigning equal coverage levels, this portion of
the analysis specifically answers the question *"whether short reads can just
as successfully as long reads, discriminate pseudogenes from their parents."*          

1x 3x 6x 10x 20x 50x        

**If possible, i'd like to add low expression values for the pseudogenes, and
much higher values for their parent genes, to see if this biases the mapping
towards the parents.** - Most likely won't have the time for this though..      

Below are the Rscripts utilised.     

```     
# Load Packages
library(polyester)
library(Biostrings)

# Load transcripts .fasta file
fasta_file <- 'custom-polyester-input-ALL-splicedtranscripts.fa'
fasta = readDNAStringSet(fasta_file)

# Assign 1, 3, 6, 10, 20, 50x read coverage to the transcripts.
readsperttx2 = round(1 * width(fasta) / 100)
readsperttx3 = round(3 * width(fasta) / 100)
readsperttx6 = round(6 * width(fasta) / 100)
readsperttx10 = round(10 * width(fasta) / 100)
readsperttx = round(20 * width(fasta) / 100)
readsperttx50 = round(50 * width(fasta) / 100)
         
# Simulate reads
simulate_experiment('custom-polyester-input-ALL-splicedtranscripts.fa',
                     reads_per_transcript=readsperttx2, num_reps=c(4,4),
                     readlen=150, fold_changes=
                     outdir='/Users/uqdbasic/Desktop/R-readsim/output-simulated-reads/1x')
```        
* Performed as above for the respective coverages.       

Didn't include an instrument specific error model, or a cDNA fragmentation
model (3' bias). Which may have a larger effect on gene-copy/pseudogene
discover than the shorter read lengths themselves. rRNA depletion effects if
sequences have rRNA homology - primer biases, pcr amplication bias, underlying sequence GC bias etc.


As illumina themselves recommened 2x75bp for their transcriptome sequencing,

```    
# let's simulate another set at 1x, 5x, 10x coverage.      
 
simulate_experiment('custom-polyester-input-ALL-splicedtranscripts.fa',
                     reads_per_transcript=readsperttx5, num_reps=c(4,4),
                     readlen=75, fold_changes=1,
                     outdir='/Users/uqdbasic/Desktop/R-readsim/output-simulated-reads/75bp-pairedend/5x')
```

### Mapping to the genome. 

I will use HISAT2 initially, generate .bam files, convert them to .bed files, then perform intersects with the novel pseudogenes in our dataset.     

Also try BWA-mem?       


#### HISAT2 

150bp paired end at 1x coverage -> map, sort, index, then convert to .bed. Intersec with novel-dRNA pseudogene set.    
* High sensitivity setting for HISAT2.        

```
HISATX=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/mm10/HISAT2-idx-grcm38/genome    

hisat -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -f -p 10 -x $HISATX -U sample_01_1.fasta,sample_01_2.fasta -S
PS-readsim-hisat-x1.sam    

samtools view -@ 10 *.sam -b -o PS-readsim-hisat-x1.bam    

samtools index *.bam   

samtools sort *.bam -o PS-readsim-hisat-x1-sorted.bam

bedtools bamtobed -i *.bam > PS-readsim-hisat-x1.bed    

# add 'chr' prefix to first bed column as HISAT index has different naming values for chromosomes (simply numbers/letter).    

less -S *.bed | sed -e 's/^/chr/' | bedtools intersect -a $PSBED -b - | awk '{print $11}' | sort | uniq | wc -l
```    
* 98% mapping rate. 78% unique mapping.    


* Remarkably, all of the transcripts were faithfully mapped and reconstructed - all 72! Even at 1x coverage.     
* This means that even short reads, with sufficient accuracy ofcourse, can
discriminate between highly similar pseudogenic transcripts, and detect their
novelty.    
* So where are they getting lost? We saw with our short-read data, that many of
the pseudogenes were present in the raw .bam files, but only 1-2 were found in
the .gtf annotation reconstructed by stringtie, likely pointing in the
direction in which they are discarded by the gene annotation programs, perhaps
because of low coverage, thus not meeting stringency cutoffs.     


#### Transcript reconstruction

Perform some transcript reconstruction with different software packages:
Stringtie, Cufflinks and Scripture.  The main drive of this section is to gain
insight into whether the transcrip assembly software is discarding the novel
transcripts - most likely due to expression levels.       

I will reconstruct the transcripts at 1x coverage, then at 20x coverage.      

Transcriptome Annotation: Gencode vm24.     


*StringTie v2.1.4*           

At 1x coverage.     

With default settings.     

```   
stringtie *.bam -G $vm24 -o ./stringtie-default.out

bedtools intersect -a $PSBED -b stringtie-default.out | awk '{print $11}' | sort | uniq | wc -l
```
* Only 24 of the transcripts have been reconstructed - not good enough.     
* Modify some parameters 

Turning off multi-mapping correction did not change anything.     

The combination of reducing -s to 1 made the largest difference: As most of the novel pseudogenes are monoexonic this is understandable. 
-s minimum reads per bp coverage to consider for single-exon transcript (default: 4.75)

Disabling multi mapping correction resulted in one additional discovery.     

Reducing -g, the maximum gap allowed between alignments to 0 also resulted in an additional assembly. 

`stringtie *.bam -G $vm24 -L -s 1 -g 0 -u | bedtools intersect -a $PSBED -b - | awk '{print $11}' | sort | uniq | wc -l`     
* 61 out of 72 transcript reconstruction with modified parameters.      
* Not sure what else I can change to stimulate the discovery of the remaining 11 transcripts?!     

At 20x coverage

* Able to reconstruction 69 transcripts at default settings 

Slightly modified settings as above leads to discovery of 71 transcripts, with
additional parameter changes not leading to the reconstruction of the
additional transcript.     
`stringtie *.bam -G $vm24 -s 1  -u  | bedtools intersect -a $PSBED -b - | awk '{print $11}' | sort | uniq | wc -l`


At 6x coverage

* Able to reconstruction 66 transcripts - so not much of an advantage with the 20x coverage.    

Slight modification as was done with the 1x set led to reconstruction of 68 of 72 transcripts.   
`stringtie *.bam -G $vm24 -s 0.5 -g 0 -u  | bedtools intersect -a $PSBED -b - | awk '{print $11}' | sort | uniq | wc -l`     


*Cufflinks version 2.2.1* - Cufflinks hasn't been updated in 7 years so keep this in mind.    

With default settings at x1 coverage. 

`
~/Bioinformatics-Programs/cufflinks-2.2.1.Linux_x86_64/cufflinks -p 10 -g $vm24 *.bam    
`      
* Default settings reconstructed 30 transcripts - so more permissive than StringTie.    


Change some paramters - reduce --min-frags-per-transfrag from 10 to 1, and enable multi-mapping correction.    
* Reconstructed 63 transcripts, slightly more than StringTie 


At 20x coverage.     

Default settings
* Reconstructs 65 transcripts - less effective than StringTie at same coverage. 

Modified settings as above
* Didn't change from 65 - same reconstruction numbers.    
* The idiosyncracies between the programs will determine why some over others
are assembled/reconstructed - and it is out of the scope at the moment to
understand exactly what these idiosyncracies are. The usual suspects, such as minimum coverage, transcript length, gap-intron size.         


#### Repeat with 75bp paired end reads - HISAT2 

* I will perform the exact same analysis as above, but this time with 75 bp
paired end reads, as this is the read length recommended for transcriptomics
work by Illumina.     


Mapping with HISAT2 - setting as used above.    
* 98.9% mapping rate, 70% mapping uniquely.        

```
hisat2 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -f -p 10 -x $HISATMOD -U
sample_01_1.fasta,sample_01_2.fasta | samtools view -@ 10 - -b | samtools sort
- -o PS-readsim-hisat-mod-75bp-x1.bam 

bedtools intersect -a $PSBED -b *.bam | awk '{print $11}' | sort | uniq | wc -l
```     
* Remarkably, all 72 transcripts mapped - indicating, even at shorter 75bp
reads, there is sufficient sequence information to discriminate between these
highly similar transcripts.     


Perform StringTie Annotations 
`stringtie *.bam -G $vm24 -s 1 -g 0 -u | bedtools intersect -a $PSBED -b - | awk '{print $11}' | sort | uniq | wc -l`
* Remarkably, at the default settings - 

