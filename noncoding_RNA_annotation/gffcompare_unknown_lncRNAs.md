# Document Outline 
Use gffcompare, alongside bedtools and bedops, to extract the known lncRNAs from the combined dRNA dataset, as well as completely unknown transcripts.

# Aims
* Have a general breakdown of the loci from which lncRNAs are being transcribed. 
* Figures for total lncRNAs in our dataset.
* Figures for all unknown transcripts within our dataset.   

# Analysis

## Pre Pseudogene,tRNA ommission
### Unknown category
#### Basic numbers
How many total 'unknown' transcripts in the dataset?
`awk '$3=="u" {print $0}' *.tmap  | wc -l`  
* 2915 total unknown 'reference' transcripts (~3219 listed with exons total)

Determine average exon density.  
`awk '$3=="u" {print $0}' *.tmap  | sort -k6 -rn | cut -f 6 | paste -sd+ | bc`
for total number of exons in the dataset (exons are stored in row 6 of the
.tmap).  
> 4542.  

Now dividing this total number by the number of transcripts. `4542/2915` = 1.6
exons per gene.  

What are the upper and lower limits?  
Highest is 14 exons (STRG.33915.1), and lowest is 1. Median is 1. **Consider
plotting this for better representation**   

A basic script to calculate min, max, avg. med. from a single column of numbers
found at
[here](https://unix.stackexchange.com/questions/13731/is-there-a-way-to-get-the-min-max-median-and-average-of-a-list-of-numbers-in)  

```
awk '$3=="u" {print $0}' *.tmap  | sort -k6 -rn | cut -f 6 | sort -n | awk
'BEGIN {c = 0; sum = 0;} $1 ~ /^(\-)?[0-9]*(\.[0-9]*)?$/ {a[c++] = $1; sum +=
$1;} END {ave = sum / c; if( (c % 2) == 1 ) { median = a[ int(c/2) ];} else
{median = ( a[c/2] + a[c/2-1] ) / 2;} OFS="\t"; print sum, c, ave, median,
a[0], a[c-1];}'
```  

Determine transcript length.   
The above script can also be used to calculate this from the 10th row in the .tmap.   
> Avg.: 2286.84	 Med.: 1747	 Min.: 200	Max: 15018.  

How many transcripts are there over 5000nt in length?
`awk '$3=="u" {print $0}' *.tmap  | sort -k10 -rn | cut -f 10 | awk '$1>5000
{print $0}' | wc -l`    
> 235 (quite uniform rise up until ~15'000nt from 5000nt).   

* I would like to bin these transcripts into ranges e.g. 1000-2000nt etc. and produce a plot.

#### Extract nucleotide sequences

Extract the nucleotide sequences of these unknown transcripts in .fasta format.  
Generate a .gtf file of the unknown transcripts
```
awk '$3=="u" {print $0}' *.tmap | cut -f 5 > ./unknown/unknown-txs.lst  
fgrep -w -f ./unknown/unknown-txs.lst gffcmp.annotated.gtf > ./unknown/unknown-txs.gtf  

gffread -w unknown-txs.fa -g $GENOME ./*.gtf -W -v 
```   
-v in gffread option (verbose/expose problems in gtf records).  
* The .tmap output lists the best reference transcripts - it does not list ALL of them.   

If we want to extract the nucleotide sequence of a specific transcript, we can
do the following.   
`seqkit grep -p "STRG.22881.11" *.fa`   STRG.9840., STRG.30833.2, STRG.30862.1


Genomic locus of the transcripts. 
* Number that are intronic, number that are intergenic.  

#### Genomic locus of origin

First, the intergenic and intronic portion of the GENCODE annotation must be
generated. Be will use bedtools for this.  

**Intergenic first**. Be sure to have the correct "genome files" for bedtools.  
`mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size
from mm10.chromInfo" > mm10.genome`     

Convert the gencode.annotation.gtf into .bed format. I have used bedops in the
past, however, it has known problems in conversion of gencode .gtf, and so
adding this additional awk script is recommended so that the following error
isn't produced. This code adds a pseudo-"transcript ID" to the .gtf when one isn't present.  
> Potentially missing gene or transcript ID from GTF attributes.   

`awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id
\"\";"; }' input.gtf | gtf2bed - > output.bed`   

Now sort with bedtools. Must include the -g option with the mm10.genome
chromsome co-ordinates in order to all correct sorting.   
`bedtools sort -g mm10.genome -i output.bed > gencode.vM24.annotation.sorted.bed`    

Complement with bedtools to isolate intergenic regions.  
`bedtools complement -i gencode.vM24.annotation.sorted.bed -g mm10.genome >
gencode.vM24.annotation.intergenic.bed`    

Another possible method for doing this, as recommended by [Dave
Tang](https://davetang.org/muse/2013/01/18/defining-genomic-regions/). However,
this is unproven yet, and it appears that the file produced by this method
ommits ChrM complements.   
`cat gencode.vM24.annotation.gtf | awk 'BEGIN{OFS="\t";} $3=="gene" {print
$1,$4-1,$5}' | bedtools sort -g mm10.genome | bedtools complement -i stdin -g
mm10.genome > testing.bed`   

Convert unknown.gtf into a .bed, in preparation for intersection.  
`cat unknown-txs.bed | bedtools sort -g
/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm24/filters/mm10.genome
> ./unknown-txs.sorted.bed` 

Intersect the two .bed files with one another to look for overlap. -wb option.
```  
INTG=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm24/filters/gencode.vM24.annotation.intergenic.bed

awk '$8=="transcript" {print $0}' ../unknown-txs.sorted.bed | bedtools intersect -a stdin -b $INTG -wb | less -S  
```    
* 2913 intergenic transcripts - quite a large number of ambigious reads. Now we must parse and parse and parse.   
* **Which genes are most proximal to each of these transcripts** -> get a list, then do some gFinder ontology investigation.  

Directionality of the transcripts.  
`awk '$8=="transcript" {print $0}' ./unknown-txs.sorted.bed | bedtools
intersect -a stdin -b $INTG -wb | cut -f 6 | sort -rn | uniq -c`  
* 57 as . (not sure what this means?).     
* 1392 on the anti-sense strand.  
* 1464 on the sense strand.  

Chromosome of origin.   
`awk '$8=="transcript" {print $0}' ./unknown-txs.sorted.bed | bedtools
intersect -a stdin -b $INTG -wb | cut -f 1 | sort -rn | uniq -c > unknown-txs_chrInfo.txt`    
* Chr4 has the most transcripts coming from it.  

#### Closest neighbouring gene   
* Closest neighbouring gene, only the first (not all occurances) - using the -t first option in bedtools.  
* Also try to output the distance to it. Might be a subtract bedtools function -> or simply subtracting with awk. Not to worry -- bedtools closest has an option for this (aaron quinlan the GOD).    

The gene_id's of the closest gene are in the 14th column, and I will merge
replicate genes using uniq -c, then cut out only the gene id's but not the
count in the adjacent column (space delimiter).     
`awk '$8=="transcript" {print $0}' ../unknown-txs.sorted.bed | bedtools
intersect -a stdin -b $INTG | bedtools closest -a stdin -b $GTF -io -t first |
cut -f 14 | uniq -c | cut -d " " -f 8 | head`    
* As the gene ontology (g:Finder) program doesn't permit the use of transcript ID's
(gene_id subset), the identifier following the gene_id must be removing (.*).
This can be done with awk and a special field delimiter - very clever way to do
it.   
`awk '$8=="transcript" {print $0}' ../unknown-txs.sorted.bed | bedtools
intersect -a stdin -b $INTG | bedtools closest -a stdin -b $GTF -io -t first |
cut -f 14 | uniq -c | cut -d " " -f 8 | awk -F '\\.*' '{print $1}' >
closest-to-unknown-intergenic.txt` 



**Intronic now.**  
Merge only the exon records/intervals from the gencode.annotation.gtf. This
will produce a single interval which is an amalgam of all exons. Everything
between these exons is intronic.    
`awk '$3=="exon" {print $0}' ../gencode.vM24.annotation.gtf | bedtools sort -g
mm10.genome | bedtools merge -i - > gencode.vM24.annotation.exon.merged.bed`

To identify intronic regions, we need to subtract the exonic regions from the
entire genic interval, which is a single, long block spanning the entire gene
body.   
`awk '$3=="gene" {print $0}' ../gencode.vM24.annotation.gtf | bedtools sort -g
mm10.genome | bedtools subtract -a stdin -b
gencode.vM24.annotation.exon.merged.bed > gencode.vM24.annotation.intronic.bed`

For a sanity check, let us intersect the intronic, and the exonic .bed files -
we should see no intersection and thus no overlap at all between the files.   
`bedtools intersect -a gencode.vM24.annotation.exon.merged.bed -b
gencode.vM24.annotation.intronic.bed`   
* No output, and so the procedure was performed properly.   

```
INTR=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm24/filters/gencode.vM24.annotation.intronic.bed

awk '$8=="transcript" {print $0}' ../unknown-txs.sorted.bed | bedtools intersect -a stdin -b $INTR -wb | less -S  
```  
* Only two transcripts originated from intronic regions - that of Snhg14, a
large snoRNA host. A complicated locus, annotated as a "fragmented locus".    
* **An interesting case study, or application of long read cDNA, dRNA nanopore
sequencing, is to find regions annotated as "fragmented loci", and attempt to
resolve some of the complications with this sequencing, along with some
short-read polishing.** 

#### Calculate GC Content 
Calculated the GC content of the sequences. Protein coding genes are known to
have an elevated GC content (due to the GC bias in codons) - this ties into the 
larger isochore structure it appears.  

`seqkit fx2tab --name --gc *.fa | sort -k 5,5 -rn > unknown-txs_GC-percentage.txt`  
* The sequence with highest GC% is 73.08%, and the lowest is 22.34%.  
* Are the ones with highest GC content, the prime candidates for potential
proteins or peptides?   

#### Gene Ontology Insights - Pathway Analysis   
Use the software g:Finder to perform the gene ontology analysis. In order to
determine whether results from the ontology analysis are relevant, I will also
extract ~2000 random gene names from our dataset, in multiple iterations, and
input them into g:Finder. The main reason for doing so, is that transcriptomic
sequencing from the mouse brain, will no doubt detect transcripts directly
involved in brain function, and so it should not be a great surprise, to find
genes that are enriched in "synaptic processes, and "neuronal system" and so
on.     


-------------------------------------------------------------------------------


## Post Pseudogene,tRNA, and ambigious scaffold ommission 

### Protein Coding Potential Set

I have realised now, after inspecting some preliminary high-confidence protein
homology targets, that many of the strongests hits are represented by
pseudogenes, which have not been included in the GENCODE annotation that I
used. As such, I must intersect the "unknown" transcripts with the
pseudogene, and tRNA annotation, then repeat the above analysis. Moreover,
inspecting some of the highest blastp and hmmer pfma hits, it appears there
are many transcripts that are mapping to the ambigious/fragmented/incomplete
scaffolds in the genome, some of which are not included in the typical
genome-chromosome files, nor within the comprehensive gencode-annotation.gtf.
As such, they are creating false positive hits, being detected as containing
pfam domains, only to then map to an already known gene when BLAT is
performed. In order to reduce any further confusion, I will omit all the the
transcripts from the unknown set which map to these scaffolds. 

#### Intersect unknown.gtf with pseudogene, and tRNA .gtf annotations.  
**Pseudogenes**   
```   
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";";
}' gencode.vM25.2wayconspseudos.gtf | gtf2bed - >
gencode.vM25.2wayconspseudos.bed    

bedtools sort -g ../vm24/filters/mm10.genome  -i
gencode.vM25.2wayconspseudos.bed > gencode.vM25.2wayconspseudos.sorted.bed   

PSBED=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm25/gencode.vM25.2wayconspseudos.sorted.bed

awk '$8=="transcript" {print $0}' ../unknown-txs.sorted.bed | bedtools intersect -a stdin -b $PSBED -wb 

awk '$8=="transcript" {print $0}' ../unknown-txs.sorted.bed | bedtools
intersect -a stdin -b $PSBED -wb | cut -f 10 | cut -d " " -f 2 | sed
's/["";]//g' > pseudogene-unknowns.lst
```   
* There are 36 transcripts which intersect with known pseudogenes. These will be removed from the unknown searches.  

**tRNAs**  
* Same procedure as above to create sorted.bed analogues.    

```
TBED=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm25/gencode.vM25.tRNAs.sorted.bed    

awk '$8=="transcript" {print $0}' ../unknown-txs.sorted.bed | bedtools intersect -a stdin -b $TBED -wb    

awk '$8=="transcript" {print $0}' ../unknown-txs.sorted.bed | bedtools
intersect -a stdin -b $TBED -wb | cut -f 10 | cut -d " " -f 2 | sed
's/["";]//g' > tRNA-unknowns.lst     
```    
* There are 155 transcripts which arise from tRNA genomic loci. This will also
be removed from the unknown searches.    

#### Remove the transcripts of pseudogenic, and tRNA origin    
`fgrep -w -v -f pseudogenised-tRNAs/pseudogene-unknowns.lst unknown-txs.gtf |
fgrep -w -v -f pseudogenised-tRNAs/tRNA-unknowns.lst - >
no-annotation-truly/true-unknown.gtf`    
* 2734 transcript remaining. We had duplicates (pseudogenised tRNAs likely).   

#### Extract .fasta nucleotide sequences   
`~/Bioinformatics-Programs/gffread/gffread -w true-unknown.fa -g $GENOME ./*.gtf -W -v`   

#### Chromosome of origin
* Found in chrmInfo.txt file.  

#### Separate out transcripts mapping to ambigious haplotigs-scaffolds 
As I peruse the output from blastp alignment, and hmmer mapping, I am noticing
that many of the transcripts with positive hits, are not originating from truly
"unknown" regions, but are instead coming from the minority scaffolds of the
genome, which are often not included in the main gencode annotation files. As
such, I have not been able to remove them from my unknown transcript set, and I
am thus creating false positives. As such, I will remove all ambiguously
mapping transcripts, and isolate them in their own set, which can be browsed
and probed further.    

* I renamed the transcripts which were previously true-unknown, to
all-ambg-unique, and left the true-unknown for the transcripts not mapping to
any ambigious scaffold, nor to pseudogenes or tRNAs.    

Create list of truly unique trascripts, with clear mapping co-ordinates. 
`cat all-ambg-unique.gtf | grep -v "chr._" | grep -v "chrUn_." | cut -f9 | cut
-d " " -f 2 | uniq | sed 's/["";]//g' > ../true-unknown.lst`    

Create the corresponding .gtf and .fa.    
```   
cat all-ambg-unique.gtf | grep -v "chr._" | grep -v "chrUn_." > ../true-unknown.gtf  
gffread -w true-unknown.fa -g $GENOME *.gtf -W -v    
```     
* 2583 transcripts are assigned strictly to high confidence chromosomes.   
* 151 transcripts were assigned to ambigious scaffolds.    

Remove the ambigious transcripts from the blastp and pfam analyses.   
```
TRU=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/true-unknown.lst

fgrep -w -f $TRU pfma.domtblout > new

cat blastp.outfmt6 | sort -k3,3 -rn | less -S
```

## Tab-delimited information file for unknown transcripts.
* I would like to create a file which contains much of the relevant information
for the category of unknown transcripts. The columns will be as following:
ID, Coverage, Number of Exons, Length, ?Closest gene?, PFAM hit, PFAM
evalue, BLAST hits, BLAST evalue, TransDecoder ORF, Trinotate detail,
InterproScan detail, and potentially, RefTSS support, In House Orthogonal
Data.     

#### Transcript level coverage 
Extract the transcript level coverage of all the unknown transcripts.   
`fgrep -w -f true-unknown.lst ../../../RC-EXT_stringtie.gtf | awk
'$3=="transcript" {print $0}' | cut -f 9 | cut -d";" -f2,3 | cut -d" " -f3,5 |
sed 's/"//g' | sort -k2,2 -rn | sed 's/;//g' > ./true-unknown-coverage.lst`    


#### Extract number of exons and transcript length
`cat gffcmp.RC-EXT_stringtie.gtf.tmap | cut -f5,6,10 | less -S | fgrep -w -f ./unknown/no-annotation-truly/true-unknown.lst - | less -S`     


#### Sorting PFAM output
Add transcript i.d. next to corresponding transdecoder ID, in an attempt to add coverage levels to a column corresponding to transcript ID.   
`awk 'F="\t" {print $1}' pfma.domtblout | sed 's/.\{3\}$//' | paste - pfma.domtblout > v2.pfma.domtblout`    

Add the respective coverage value for each transcript within the PFAM output. Of note, one must ensure that each duplicate transcript ID, which may have a unique transdecoder .p ID, also recieves the same coverage value adjacent to it. 
* This is what awk was created for, the forums say, [here]( https://askubuntu.com/questions/890557/awk-compare-2-files-and-print-columns-from-both-files), and [here](https://unix.stackexchange.com/questions/134829/compare-two-columns-of-different-files-and-print-if-it-matches).  
* Using gawk, we can read read the values in the first column of the first
  file, scan the strings in the column of the second file, and if we recieve a
  match, to print out the second column of the first file, the one containing
  the coverage values, next to the correct transcript ID. This was a
  challenging one!   

`gawk 'NR==FNR {a[$1][$2]++; next} $1 in a {for (x in a[$1]) print x, $0}' OFS="\t" $lst v2.pfma.domtblout | sort -k9,9 -g |`   

Extract only those with an e-value score over 1e-3.   
`awk '$9 < 1e-03 {print $0}' v3.pfma.dombtblout > v3.threshold.pfma.dombtblout`   
* There are 134 unique PFAM hits over threshold.   

#### Sorting BLASTP output    
Add transcript i.d. next to corresponding transdecoder ID.   
`awk 'F="\t" {print $1}' blastp.outfmt6 | sed 's/.\{3\}$//' | paste - blastp.outfmt6 > v2.blastp.outfmt6`   

Add the respective coverage as was performed above.  
```
lst=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/true-unknown-coverage.lst  

gawk 'NR==FNR {a[$1][$2]++; next} $1 in a {for (x in a[$1]) print x, $0}'
OFS="\t" $lst v2.blastp.outfmt6 > v3.blastp.outfmt6   
```   

Filter to remove all hits with e-values lower tham 1e-3.   
`awk 'F="\t" $11 < 1e-03 {print $0}' v3.blastp.outfmt6 | sort -k13,13 -g > v3.threshold.blastp.outfmt6`    

#### Add coverage information to transdecoder.bed   
`FA=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/true-unknown.fa`   

`gawk 'NR==FNR {a[$1][$2]++; next} $1 in a {for (x in a[$1]) print x, $0}'
OFS="\t" $lst true-unknown.fa.transdecoder.bed | sort -k1,1 -rn >
v3.true-unknown.fa.transdecoder.bed`   

#### Add coverage information to signalp, and tmhmm output.   
```   
awk 'F="\t" {print $1}' signalp.out | sed 's/.\{3\}$//' | paste - signalp.out > v2.signalp.out   

gawk 'NR==FNR {a[$1][$2]++; next} $1 in a {for (x in a[$1]) print x, $0}'
OFS="\t" $lst v2.signalp.out | sort -k1,1 -rn > v3.signalp.out   

awk 'F="\t" {print $1}' tmhmm.out | sed 's/.\{3\}$//' | paste - tmhmm.out > v2.tmhmm.out   

gawk 'NR==FNR {a[$1][$2]++; next} $1 in a {for (x in a[$1]) print x, $0}'
OFS="\t" $lst v2.tmhmm.out | sort -k1,1 -rn > v3.tmhmm.out   
```

#### Create a 'super' information file containing all the relevant information for each transcript. 
* As the ORF/Protein prediction works on the transdecoder ID's with the *.p2*
suffix, I have to use this information as the main identifier when collating
this tab-separated file.   
* The merged/collapse 'super' file will only contain blastp + pfam hits over
threshold of e-03, and only transmembrane predictions with a positive hit.     
* Naming each file "pre*suffix* before the merge.    
* **One HUGE caveat, I must discuss. As TransDecoder breaks each transcript in
separate subsequences based on it's gene structure e.g. 3'UTR, exon, CDS etc.,
the transcript as given sub-ID's such as STRG101.p1, .p2, .p3 and so on.
Running these 'peptides' through blastp and pfam, can result in multiple
positive hits for each transcripts, as the transcript could be a multi-domain
protein, with repeating domains all clustered together. As such, it is
impossible to create a 1:1 correspondance between all of the files, as one of
the programs will inevitibly score/predict regions which the other(s) will not,
leading to discordance between all the transcript sub-IDs. e.g. PFAM may have
10 hits for STRG101.p1, while BLASTP, only has 2. Because of this difference,
exact matches between the sub-sequences is not possible. HOWEVER, the following
programs, will output ANY of the matches between programs, and as such, if one
requires more specific information for the transcript at hand, one should go
directly to the respective output files. The 'super-file' is useful as a basic
reference which can quickly tell you if a transcript has multiple matches
between programs/databases, without having to constantly bounce back and forth
between the individual output files.**     

Coverage : Length : Exons : Transcript-ID : TransdecoderID : PFAM hits : BlastP Uniref90 Hit : SignalP : Transdecoder ORF : Transmembrane Domain : 
OrthoDB : DFAM : Poly-A cluster-site


**Testing, hacking**  
`awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1 in a)?a[$1]:"NA"}' pfma.domtblout blastp.outfmt6 ` 

Turn whitespaces in a file into tab - essentially delimiting different columns of data.. 
`cat tmp | awk -v OFS="\t" '$1=$1' | less -S` ***

The progressive intersects work - as below 
`awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1 in a)?a[$1]:"NA"}'
../pre-tmhmm tmp | awk -v OFS="\t" '$1=$1' | less -S`    

1) Create a list with all Transdecoder IDs in the first column.  
`cut -f 2 gene_trans_map > Transcoder-IDs.tab`   

2) Parse the Transdecoder file to extract only the relevant pieces of
information from it - Transdecoder ID, ORF-type, length, strand, and score. The
following command will progressively cut columns and change delimeters to create
the needed format.     
`cat v3.true-unknown.fa.transdecoder.bed | cut -f5 | cut -d";" -f 2,3 | cut -d,
-f1,2 | sed 's/.*~//' | sed 's/;/\t/g' | sed 's/,/\t/g' > `

3) Intersect/merge all the relevant files to create the final 'super file'.    
`awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1 in a)?a[$1]:"NA"}'
pre-pfma Transcoder-IDs.tab | awk 'NR==FNR {a[substr($1,1)]=$0; next} {print
$0,($1 in a)?a[$1]:"NA"}' pre-blastp - | awk 'NR==FNR {a[substr($1,1)]=$0;
next} {print $0,($1 in a)?a[$1]:"NA"}' ../pre-tmhmm - | awk 'NR==FNR
{a[substr($1,1)]=$0; next} {print $0,($1 in a)?a[$1]:"NA"}'
Trandecoder-bed-chopped.tab - | `   
**A monstrosity of a command**    (This is why I need to learn more python and SQL!!!!).    

Signal-P information? 

4) Add the coverage levels of the transcripts to the file.     
* Performed as outlined above.    


#### Determine whether the high confidence transcript have a poly-A cluster site in them.     
* Could I aslo extend the read *x* bases in each direction or no?    
* Create .bed file from transcripts. There are 794 transcripts, so there should be 794 entries in the bed file.     

`awk '$3=="transcript" {print $0}' true-unknown.gtf | grep -f $TMP -w - | gtf2bed - > Superfile-IDs.bed`    

Sort bed file.  
```  
GENOME=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm24/filters/mm10.genome   
bedtools sort -g $GENOME -i Superfile-IDs.bed > Superfile-IDs_sorted.bed     
```   

Append "chr" to the beginning of the poly-A.bed, then intersect with Superfile.bed 
```   
POLYA=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/Poly-A-Sites/atlas.clusters.2.0.GRCm38.96.bed    

awk '$1="chr"$0' $POLYA    

awk '$1="chr"$0' $POLYA | bedtools intersect -a Superfile-IDs_sorted.bed -b stdin -wb   
```    
* 422 transcripts have PolyA cluster support. Interesting news no doubt.   

I will use this file to append the Superfile.tab with Poly-A cluster details.  
`awk '$1="chr"$0' $POLYA | bedtools intersect -a Superfile-IDs_sorted.bed -b stdin -wb | cut -f 10,12-21 | cut 
-d";" -f 1,6 | sed 's/transcript_id "//g' | cut -f1,11 | sed 's/;//' | sed 's/"//' > Superfile-PolyA-pre.tab`     

`cat $POLS | awk '{print $2,$1}' | awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1 in a)?a[$1]:"NA"}' - Superfile_Cov.tab ` **work in progress**

Append Superfile with Poly-A information.   
`awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1 in a)?a[$1]:"NA"}' $POLS tmp-Super.tab > tmp2.Super.tab`    


#### Add SignalP information
* Performed as described above for tmhmm data, with the inclusion of `grep -v "#"` to remove hash characters in the input file.   

`gawk 'NR==FNR {a[$1][$2]++; next} $1 in a {for (x in a[$1]) print x, $0}' OFS="\t" $lst tmp3.Super.tab | sort -k1,1 -rn > ./TransDecoder/Superfile_Cov.tab`    

#### Search for orthologues using OrthoDB database, and BLAST search.    
* Should I searchq against OrthoDB with the transdecoder "peptide", or with the original FASTA sequence? Both? 
* I went with the peptide-orfs predicted by transdecoder, as these were the
same ones used as input for the pfam, and uniref90 searches.   

Create OrthoDB blast database.  
`makeblastdb -in $FASTA -input_type fasta -dbtype prot -title orthodb10v1
-parse_seqids -out
/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Proteomes/OrthoDB/`

```   
ORTHO=/media/labpc/Disk-2/Bioinformatics-Computational/blastdb/OrthoDB/OrthoDB

PEP=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/TransDecoder/true-unknown.fa.transdecoder.pep

blastp -db $ORTHO -query $PEP -max_target_seqs 1 -outfmt 6 -evalue 1e-3
-num_threads 8 >
/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/orthodb.outfmt 
```     
* 510 have detectable orthology to existing proteins - hmmmm.   

Add Output information to SuperFile.   
```   
OT=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/TransDecoder/orthodb.outfmt   

awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1 in a)?a[$1]:"NA"}'
./pre-signalp ./TransDecoder/Superfile.tab | awk 'NR==FNR {a[substr($1,1)]=$0;
next} {print $0,($1 in a)?a[$1]:"NA"}' $OT - > tmpSUP
```  



#### Run hmmscan against the DFAM database     
* Downloaded the dfam.hmm and associated .hmm files from
[here](https://www.dfam.org/releases/Dfam_3.3/infrastructure/dfamscan/).    

```    
DHMM=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/DFAM/hmm/Dfam.hmm    

hmmsearch --cpu 2 --domtblout dfam.domtblout $DHMM $CDS   
```    
* Far too many matches, as I searched the entire list of *best-orfs* from the
TransDecoder output. As such, many of these ORFs are not deemed to be
significant in any of the analyses (blastp, pfam etc.). 
* Instead, I will extract all of the IDs from the superfile which have atleast
one substantial hit - and then rescan these against DFAM. This can be done with seqtk.       

``` 
seqtk subseq *.cds tmp.lst > Superfile-best-Orfs.fa   

CDS=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/TransDecoder/Superfile-best-Orfs.fa    
q
hmmsearch --cpu 8 --domtblout dfam.domtblout $DHMM $PEP   
```     

Filter results to remove any hits with an e-value less than 1e-3.   
`awk '$7 <= 1e-3 {print $0}' dfam.domtblout > v2.threshold.dfam.domtblout`    

Append to SuperFile. Some of the columns with additional information were omitted, for brevity sake.    
`less -S ./TransDecoder/v2.threshold.dfam.domtblout | sort -k1,1 | awk '{print
$1,$2,$3,$4,$5,$6,$7}' | awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1
in a)?a[$1]:"NA"}' - tmpSUP | less -S`    

#### Add finishing touches to SuperFile now that all the required information has been added to it    
`awk 'F="\t" {print $1}' tmp2SUP | sed 's/.\{3\}$//' | paste - tmp2SUP | gawk
'NR==FNR {a[$1][$2]++; next} $1 in a {for (x in a[$1]) print x, $0}' OFS="\t"
$lst - > ./TransDecoder/Superfile_Cov_final.tab`   

#### Number of exons, and transcript length 
`awk '{print $5, $6, $10}' $tmap | fgrep -w -f tmp-Stringtie-Ids.lst - | awk '{print $2}' | sort | uniq -c | head -n 20`   
* Majority are monoexonic (579), 114 have 2, 58 have 3, 24 have 4, 10 have 5, 4 have 6, and 2 have 7.  
* My questions is, of the known protein coding genes, how many are monoexonic?
And so what is the underlying "expection" value or distribution. How does that
compare to this dataset?       


* awk $2="" represents the column that should NOT be printed.     

`awk '{print $2,$0}' Superfile_Cov_final.tab | awk 'NR==FNR {a[$1][$2]++; next}
$1 in a {for (x in a[$1]) print x, $0}' OFS="\t" ../Superfile_exons-length.lst
- | awk '{$2=""; print $0}' - | awk '{print $3,$0}' - | awk 'NR==FNR
{a[$1][$2]++; next} $1 in a {for (x in a[$1]) print x, $0}' OFS="\t"
Superfile_exons-length2.lst - | awk '{$2=""; print $0}' - >
../Superfile_Cov_finale.tab`     



#### Closest neighbouring gene, and distance to it 
* May be potentially lengthened 5-3' UTRs, or even uORFs. 
* My curiousity exists around the chance that a lengthened UTR is in fact
protein coding, or contains a viable ORF - my suspicion is that this is
unlike.y, as the UTRs may instead carry varying regulatory signals rather than
bona-fide functional domains...  



#### Run BLASTN against RefSeq known genes    
* I will also run this against the candidate non-coding dataset.    
* [http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/](http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/) 
* I am only using the known, curated genes (NM* prefix) for I am searching
against this database in order to discover primarily pseudogenes, but also
potentially highly similar gene copies.     

`~/Bioinformatics-Programs/gffread/gffread -w mm10.knownGene.fa -g ../../mm10/mm10.fa ./mm10.knownGene.gtf -W -v`    

The following code may be able to assist in replacing the "transcript_id" with the respective "gene_id" from the annotation.gtf .     
`awk -v FS="\t" 'NR==FNR {match($9, /transcript_id "([^"]+)"/ , t); match($9,
/gene_name "([^"]+)"/, g); transcript[t[1]]=g[1]; next} {match($0, />([^ ]+)/,
t); gsub(t[1], transcript[t[1]]); print}' transcripts.gtf input.fasta >
output.fasta`     

Build blast database.    

```   
REFGENE=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/RefSeq/mm10/mm10.knownGene.fa     

makeblastdb -in $REFGENE -input_type fasta -dbtype nuc -title refseqknowngene
-parse_seqids -out /media/labpc/Disk-2/Bioinformatics-Computational/blastdb/refseq/knownGene/knownGene/   
```     

Run blastn against the primary nucleotide sequence of the ORFs predicted by transdecoder. 
```
REFSEQ=/media/labpc/Disk-2/Bioinformatics-Computational/blastdb/refseq/knownGene/knownGene     

lst=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/TransDecoder/tmp2-Stringtie-Ids.lst

seqkit grep -f $lst $FA > ./blastn/true-unknown-pc-candidates-fullfasta.fa    

blastn -db $REFSEQ -query ./*.fa -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > true-unknown.fa.transdecoder-blastn.outfmt
```     
* 986 potential hits spread across - 479 unique transcripts. 
* As I am looking for bona-fide paralogues, and pseudogenes, I am going to
prioritise, and thus only focus on, the hits with maximal alignmnent length,
which span a large fraction of the query. By doing so, I will eliminate weak
partial matches, and hits to tranposable elements.     

* Also, in order to increase specificity, I will also intersect the unknown
candidates, with the refseq annotation, as I have observed that some
transcripts are already known.

#### Remove RefSeq annotations and TE associated proteins from SuperFile    

Remove refseq hits.    
* Use the suffix ${f}-2. to designate the newly filtered files.    
* Use a variation of the below on every file within the folder, then move old files to ./old-superfiles/pre-refseq-removal/     

`cat ./RefSeq-intersect/*.lst | fgrep -w -v -f - tmp2-Stringtie-Ids.lst > tmp2-Stringtie-Ids-2.lst`

Remove TE hits.     
* Manually parse v3.threshold.pfma.dombtblout PFAM entries to determine which ones are from transposable element domains.   
* Output these PFAM IDs into their own tab-separated file, to be used for intersecting against the main Superfiles.    
* 8 entries - .23963.1 was removed as it overlapped existing RefSeq annotation.     
 
RVT_1   PF00078 PF00078.29      STRG.10477.1
Transposase_22  PF02994 PF02994.16      STRG.4136.1
Gag_p30 PF02093 PF02093.18      STRG.23963.1
Gag_p24 PF00607 PF00607.22      STRG.2467.1
Tnp_22_dsRBD    PF17490 PF17490.4       STRG.4136.1
TLV_coat        PF00429 PF00429.21      STRG.18407.1
Transposase_22  PF02994 PF02994.16      STRG.23600.1
Gag_p24_C       PF19317 PF19317.1       STRG.20762.1


Several iterations of the below - all files have the ${f}-3 suffic attached to them to signify this removal.    
`cat ./pfam/TE-domains.tab | cut -f4 | fgrep -w -v -f - gene_trans_map-2 > gene_trans_map-3`     


### Non-Coding Set 

#### Closest neighbouring gene, distance to it, number of exons, transcript length      
* This will aid in identifying which non-coding transcripts may be functions as
enhancers, and/or in *cis*. Extensions of UTRs of already annotated elements,
or themselves, arising from promoter regions (if overlapping a H3K4me3 mark?)    
* Is there a relationship between exon number and coverage/count?     

fgrep -v -w -f with the Superfile-IDs against the main .lst file in order to
extract the transcript IDs of the non-coding candidates.   
* 1792 non-coding candidates.     

Create a .gtf file in order to generate .bed file for bedtools intersect.   
`fgrep -w -f true-unknown-noncoding.lst ../true-unknown.gtf > true-unknown-noncoding.gtf`   

Convert to .bed.     
`cat true-unknown-noncoding.gtf| gtf2bed - > true-unknown-noncoding.bed`   

Sort with bedtools.   

```
BEDGENOME=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm24/filters/mm10.genome   

bedtools sort -g $BEDGENOME -i true-unknown-noncoding.bed > true-unknown-noncoding.sorted.bed   
```    

Bedtools closest. Use -D ref option to input distance with respect to the
genome, but may need to use option -D a to maintain strand information. -io
option to ignore overlaps, and -t to list the first closest element.  
```   
GENBED=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm24/filters/gencode.vM24.annotation.sorted.bed    

cat $GENBED | less -S | awk '$8 == "gene" {print $0}' > tmp.bed

awk '$8 == "transcript" {print $0}' true-unknown-noncoding.sorted.bed |
bedtools closest -a stdin -b tmp.bed -io -t first -D ref |   
```    

Parse the output in order to extract the transcript ID of the query, along with
the type of gene, the gene name, and the distance to, the closest gene.     
`awk '$8 == "transcript" {print $0}' true-unknown-noncoding.sorted.bed |
bedtools closest -a stdin -b tmp.gtf -io -t first -D ref | cut -f 10,6,20,21 |
awk  '{print $3,$15,$17,$NF}' | sed 's/"//g' | sed 's/;/\t/g' >
closest-to-true-unknown-noncoding.tab`    
* wc -l is still 1792 (good).      
* Use sort -Vd for absolute sorting of numbers.    
* *Entries which are not associated with a "gene" will be run through samtools
once again, with the unmodified gencode.gtf annotation - I included only
entries which were listed as "genes" initially.*   

**A slight modification is now required as we need the transcripts strand, including the closests elements strand**

`awk '$8 == "transcript" {print $0}' true-unknown-noncoding.sorted.bed |
bedtools closest -a stdin -b tmp.gtf -io -t first -D ref -g
/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm24/filters/mm10.genome
| awk -F "\t" '{print $10,$6,$20,$21,$16}' - | awk  '{print
$2,$11,$15,$17,$(NF-1), $NF}' | sed 's/"//g' | sed 's/;/\t/g' | less -S` 


#### How many are within 5kb(-.+) of another gene?   
* This could point in the direction of extended UTRs, and perhaps enhancers (if methylation signal is present).    
* Use this special awk feature to recognize absolute numerical values.      
`less -S closest-to-true-unknown-noncoding.tab | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($4) < 5000) print $0}' | wc -l`    
* 924     

Of these, which category of gene are they associated with?    
> 134 unnamed (what could this be), 27 TEC, 1 TR_V, 95 lncRNA, 9 miRNA, 2 misc_RNA, 79 processed
pseudogene, 541 protein coding, 1 pseudogene, 1 rRNA, 7 snRNA, 5 snoRNA, 2
transcribed-processed pseudogene, 2 transcribed-unitary pseudogene, 5
transcribed-unprocessed pseudogene and 12 unprocessed pseudogene. **This has been updated now that corrections have been made**.       

Length distribution to closest gene.   
Average:
`less -S closest-to-true-unknown-noncoding.tab | awk -F'\t' 'function
abs(x){return ((x < 0.0) ? -x : x)} {if (abs($4) < 5000) print $0}' | cut -f4 |
awk 'function abs(x){return ((x < 0.0) ? -x : x)} {sum+=abs($1)} END { print
"Average = ",sum/NR}'`   
* 1608.64 

#### Determine whether the -1 distance unamed transcripts are closest to "UTR" category in the gencode annotation 
* There are 134 unnamed transcripts, all with are a distance of one bp
downstream a genetic element - my suspicion is that the 'UTR' category in the
gencode annotation is not encompassed by the 'transcript' category which I used
during the bedtools closest operation - as such, these transcripts are not
being assigned to any 'transcript' annotation. The other alternative is that
the parsing step did not extract the associated details from the closest
intersect. 

```
GENGTF=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm24/gencode.vM24.annotation.gtf   

awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' $GENGTF | gtf2bed - > tmp-gencode-UTR.bed   

cat tmp-gencode-UTR.bed | awk '$8=="UTR" {print $0}' > tmp-gencode-UTR.bed    

awk '$8 == "transcript" {print $0}' true-unknown-noncoding.sorted.bed |
bedtools closest -a stdin -b tmp-gencode-UTR.sorted.bed -io -t first -D ref |
less -S     

awk '$8 == "transcript" {print $0}' true-unknown-noncoding.sorted.bed |
bedtools closest -a stdin -b tmp-gencode-UTR.sorted.bed -io -t first -D ref |
less -S | fgrep -w -f tmp.lst - | wc -l
```

* It would appear, that the respective transcripts do not in fact, have any
existing annotations assigned to them, for whatever reason, bedtools closest is
not reporting any additional details, aside from the fact that they are -1 away
from another genetic element - suprisingly also, all such transcripts originate
from chr 18 and 19.    

#### Number of exons 

**Pasted from above section** 

Determine average exon density.  
`awk '$3=="u" {print $0}' *.tmap  | sort -k6 -rn | cut -f 6 | paste -sd+ | bc`
for total number of exons in the dataset (exons are stored in row 6 of the
.tmap).  
 

tmap=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-read/RC-and-EXT/gffcompare/gffcmp.RC-EXT_stringtie.gtf.tmap

awk '{print $5, $6, $10}' $tmap | fgrep -w -f true-unknown-noncoding.lst -  | wc -l  
* Sanity checked - 1792 transcripts (original)  

Exon density among transcripts.   
`awk '{print $5, $6, $10}' $tmap | fgrep -w -f true-unknown-noncoding.lst -  | awk '{print $2}' | sort | uniq -c | head -n 20`
* The overwhelming majority are monoexonic - 1406. The rest:1 transcript has 11 exons, 240 have 2, 97 have 3, 30 have 4, 15 have 5, and 3 have 6.    
* An intersection between monoexonic, low expression, and close proximity to another gene, will be an interesting take.   


#### Poly-A Cluster support 
```    
POLYA=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/Poly-A-Sites/atlas.clusters.2.0.GRCm38.96.bed     

awk '$1="chr"$0' $POLYA | bedtools intersect -a tmp.bed -b stdin -wb | cut -f
10,12-21 | cut -d";" -f 1,6 | sed 's/transcript_id "//g' | cut -f1 | sort |
uniq | wc -l    

awk '$1="chr"$0' $POLYA | bedtools intersect -a tmp.bed -b stdin -wb | cut -f 10,12-21 | cut -d";" -f 1,6 | sed 's/transcript_id "//' | cut -f1,11 | s
ed 's/;//' | sed 's/"//' > true-unknown-noncoding-polyA-support.lst    
```    
* One has to keep in mind, that many of these have multiple poly-A clusters
within them - utilising different poly-adenylation sequences. If I am going to
dig deeper into the transcript, I will have to go back into the primary files.   
* 788 have poly-A support - of the 1792. This is very interesting to see.    

#### TSS database support    
* Merge refTSS, and EPD Non-coding, Protein-coding datasets.  

```   
sort -k1,1 -k2,2n on each .bed file   
cat tmp.reftss tmp.nc tmp.pc | bedtools merge -i - > EPD-REFTSS-merge.bed   
```    
* 170209 transcription start site entries.    

I will extend each transcript by 500nt on either side if it is more than 500nt
away from the closest neighbouring element element. This should increase the
chances of covering a bona-fide transcription start site, and also provide
hints as to whether the transcripts are being omitted as byproduct of Pol-II
binding and thus promoter activity.      

How many transcripts are within 500nt of another genetic element?   
`less -S tmp-new-superfile.lst | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($5) < 500) print $0}' - | less -S | wc -l`    
* 251. Of these, 199 are protein coding, 19 are lncRNAs, 3 are miRNAs, 3 are TEC, 1 is snRNA, 25 are pseudogenic.   

**INCOMPLETE**


#### Compile another 'Superfile' for the non-coding dataset 
Information in columns is as follows:
*Stringtie ID : Coverage : Closest-gene type : Closest-gene name : Closest-distance to : Exon count : Transcript length: Poly-A cluster support.*    
* Superfile for entires that do not have a closest annotation will be the same as above, but without the "Closest-" information.     
```   
awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1 in a)?a[$1]:"NA"}'
closest-to-true-unknown-noncoding.tab true-unknown-noncoding-coverage.lst | awk
'{print $1,$2,$4,$5,$6}'    
```    
* I will be removing the Stringtie ID from the columns when it is added
uneccessarily, hence why I have omitted a column with awk when printing the
output.   

Exclude the entries that do not have a "closest" annotation ascribed to them. 
```
awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1 in a)?a[$1]:"NA"}'
closest-to-true-unknown-noncoding.tab true-unknown-noncoding-coverage.lst | awk
'{print $1,$2,$4,$5,$6}' | awk '$3 == "-1" {print $0}' > tmp_-1_.lst`    

awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1 in a)?a[$1]:"NA"}'
closest-to-true-unknown-noncoding.tab true-unknown-noncoding-coverage.lst | awk
'{print $1,$2,$4,$5,$6}' | awk '$3 != "-1" {print $0}' > tmp-superfile.lst   
```    

Process the remainder 

Entries with closest annotation   
`awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1 in a)?a[$1]:"NA"}'
true-unknown-noncoding-exon-number-length.lst tmp-superfile.lst | awk '{print
$1,$2,$3,$4,$5,$7,$8}' | awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1
in a)?a[$1]:"NA"}' true-unknown-noncoding-polyA-support.lst - >
tmp-superfile-2.lst`     

**Entries without** 
`awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1 in
a)?a[$1]:"NA"}' true-unknown-noncoding-exon-number-length.lst tmp_-1_.lst | awk
'{print $1, $2, $3, $5, $6}' | awk 'NR==FNR {a[substr($1,1)]=$0; next} {print
$0,($1 in a)?a[$1]:"NA"}' true-unknown-noncoding-polyA-support.lst - >
tmp_-1_2.lst`     
* 53 of the transcripts have Poly-A cluster support. Quite promising and interesting.     

**Slight Detour**
* Not being content with the -1 bedtools entries, which supposedely have no
genome annotations assigned to them, yet are shown to be only -1 away from an
unnamed genetic element, I wanted to further investigate why this was
happening, and perhaps if it could be fixed. 
* My theory was that the bedtools closest did not execute correctly, perhaps
due to the absence of a -g .genome file, which upon inclusion, appeared to
fix the problem, resulting in the expected outputs. 
* Below I have redone the closest, with the transcripts which were previously
assigned to a -1 category.   

```   
GENBED=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm24/filters/gencode.vM24.annotation.sorted.bed   
awk '$8 == "gene" {print $0}' $GENBED > tmp-gen.bed   

less -S closest-to-true-unknown-noncoding.tab | awk '$2 == "-1" {print $0}' |
cut -f1 | fgrep -w -f - true-unknown-noncoding.sorted.bed | awk '$8 ==
"transcript" {print $0}' |  bedtools closest -a stdin -b tmp-gen.bed -io -t
first -D ref -g
/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm24/filters/mm10.genome
| cut -f 10,6,20,21 | awk  '{print $3,$15,$17,$NF}' | sed 's/"//g' | sed
's/;/\t/g' | less -S
```
* All 134 entries now show the expected output, with a gene name, gene type, and distance to it, in the respective columns.     
* I shall append these data to the Superfile once the additional information has been added to it.    

Append existing Superfile with additional data from above. 
`awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1 in a)?a[$1]:"NA"}'
tmp_-1_2.lst true-unknown-noncoding-coverage.lst | awk '$3 != "NA" {print $0}'
| awk '{print $1,$2,$4,$5,$6}' | awk 'NR==FNR {a[substr($1,1)]=$0; next}
{print $0,($1 in a)?a[$1]:"NA"}' true-unknown-noncoding-exon-number-length.lst
- | awk '{print $1,$2,$3,$4,$5,$7,$8}' | awk 'NR==FNR {a[substr($1,1)]=$0;
next} {print $0,($1 in a)?a[$1]:"NA"}' true-unknown-noncoding-polyA-support.lst
- | less -S | cat - tmp2-superfile.lst > tmp-new-superfile.lst`  

#### Sorting out strandedness for non-coding transcript candidates 
* One of the greatest issues with the above "Super" information files, is that
they do not contain any information on both the transcript candidates strand,
nor the strand of it's closest genetic element. This presents an issue,
particular when one is looking to break the transcripts down into more
meaningful sub-categories of potential elements, such as enhancer, UTR
extensions, promoter associated transcripts (promoter upstream
transcripts/PROMPTS), for one requires the strand information, if one is
looking to infer whether something is acting up or downstream, and so on.    
* I have no sought to add this valuable information into the super-file. 

`awk 'NR==FNR {a[substr($1,1)]=$0; next} {print $0,($1 in a)?a[$1]:"NA"}'
true-unknown-noncoding-exon-number-length.lst
closest-to-true-unknown-noncoding-stranded.tab | awk '{print
$1,$2,$3,$4,$5,$6,$8,$9}' - | awk 'NR==FNR {a[substr($1,1)]=$0; next} {print
$0,($1 in a)?a[$1]:"NA"}' - true-unknown-noncoding-coverage.lst | awk '{print
$1,$2,$9,$10,$4,$5,$6,$8,$7}' - | awk 'NR==FNR {a[substr($1,1)]=$0; next}
{print $0,($1 in a)?a[$1]:"NA"}' true-unknown-noncoding-polyA-support.lst - |
less -S`     
* This larger file can now be (more) meaningfully parsed into the appropriate sub-categories.     



#### RFAM search 

Following the tutorial/guide as described [here](https://docs.rfam.org/en/latest/genome-annotation.html).    
Download current Rfam.cm and Rfam.clanin from ftp server.    

Compress Rfam.cm with Infernal.    
`cmpress /media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/RFAM/Rfam.cm`    

Determine how many nucleotides ("genome size") one is searching against. 
```   
NCFA=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/non-coding-set/true-unknown-noncoding.fa   

CLAN=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/RFAM/Rfam.clanin     

CM=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/RFAM/Rfam.cm    

esl-seqstat $NCFA    
``` 
> Format:              FASTA
Alphabet type:       DNA
Number of sequences: 1792
Total # residues:    3276427
Smallest:            200
Largest:             12100
Average length:      1828.4  

* 3276427 residues = "genome size"       
* To find the "Z" score, we need to multiply the number of residues by 2 (as we're searching both strands), then divide by 1'000'000 = 6.552854 

Run cmscam.   
`cmscan -Z 6.552854 --cut_ga --rfam --nohmmonly --tblout nc-rfam.tblout --fmt 2 --clanin $CLAN $CM $NCFA > nc-rfma.cmscan`     
* One result - mir-684 on a potential pseudogene.     

Testing on all of the "unknown.fa" (including potentially protein coding).     
```    
PCFA=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/true-unknown.fa        

cmscan -Z 11.413686 --cut_ga --rfam --nohmmonly --tblout pc-rfam.tblout --fmt 2 --clanin $CLAN $CM $PCFA > pc-rfma.cmscan
```     
* A few hits, some potentially interesting LSR-rRNA hits, and perhaps a few miRNA hits.    
* Do some additional, small RNA specific scans.    


#### Layer GENCODE cis-regulatory element tracks over non-coding transcripts 
* In house datasets which may be relevant include: H3K27ac (enhancer), p300 & CPB chip.     
Download CRE.bigbed file. 
`wget http://hgdownload.soe.ucsc.edu/gbdb/mm10/encode3/ccre/encodeCcreCombined.bb`    

Convert to a .bed file format with UCSC utilities "bigBedToBed".    
`~/Bioinformatics-Programs/UCSC-utils/bigBedToBed ./encodeCcreCombined.bb  encodeCcreCombined.bed`     

Split into enhancer, promoter, and CTCF-DNASE subfiles. More information can be
found
[here](https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1132306051_MHAcHICZlZKVnC2OBI6W0IBAnItQ&db=mm10&c=chr12&g=encodeCcreCombined).           
```   
less -S encodeCcreCombined.bed | awk '$10 ~ ".ELS" {print $0}' - > encodeCcreCombined-ELS.bed    
less -S encodeCcreCombined.bed | awk '$10 ~ "PLS" {print $0}' - > encodeCcreCombined-PLS.bed    
less -S encodeCcreCombined.bed | awk '$10 ~ "DNase." {print $0}' - > encodeCcreCombined-DNASE-H3K4me3.bed    
less -S encodeCcreCombined.bed | awk '$10 ~ "CTCF-only" {print $0}' - > encodeCcreCombined-CTCF-only.bed     
``` 

Intersect ENCODE-CRE.bed with noncoding candidates .bed. Ensure that the ENCODE.bed is sorted correctly.    
```    
CREBED=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/ENCODE-CRE/encodeCcreCombined.bed    
NCBED=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/non-coding-set/true-unknown-noncoding.sorted.bed

awk '$8 == "transcript" {print $0}' $NCBED | bedtools intersect -a stdin -b local.bed -wb     
awk '$8 == "transcript" {print $0}' $NCBED | bedtools intersect -a stdin -b local.bed -wb | less -S | cut -f10 | cut -d";" -f1 | sort | uniq | wc -l   
```    
* 699 transcripts have at least one sort of cis-regulatory element mark associated with them. 
* The challenge now is to intelligently, and carefully, divide the transcripts
into discrete categories based on additional sequence features, and then
re-intersect these subsets with particular CRE marks such as promoter,
enhancer, and so on.    

General overview of the types of cis-regulatory marks that are intersecting.   
`awk '$8 == "transcript" {print $0}' $NCBED | bedtools intersect -a stdin -b
local.bed -wb | less -S | cut -f10,20 | cut -d ";" -f 1,6 | cut -f2 | sort |
uniq -c`    
* Many ELS (mostly distal), with far fewer promoter associated marks.   

#### Transcripts which DO NOT have CRE signature support
* All the transcripts which don't overlap CRE signatures. 
* These will be parsed into potential enhancers (proximal and distal), prompts, 3' run ons, and bonafide intergenic transcripts.     

`awk '$8 == "transcript" {print $0}' $NCBED | bedtools intersect -v -a stdin -b
local.bed -wb | less -S | cut -f10 | cut -d";" -f1 | sort | uniq | awk '{print
$2}' | sed 's/"//g' | fgrep -w -f - ../tmp-newest-superfile.lst >
../NO-CRE-support/non-coding-NO-CRE-superfile.tab`    
* 1093 transcripts in this set - some will consist of pseudogenes and will thus be parsed out later on.      
* Quite a subsantial amount more here than in the CRE overlap set.      


#### RefSeq annotation overlap with non-coding transcripts without CRE signatures    

For completeness, remove transcripts which overlap with RefSeq annotation
(knownGene version?) - as I am finding many high confidence transcripts have
already been annotated, but not added to the GENCODE vm24 version I have been
using.  

Create a basic tab separated file which lists the transcript id, alongisde the Gene-ID from RefSeq which it overlaps with.     
```   
REFGTF=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/RefSeq/mm10/mm10.refGene.gtf     

less -S ../NO-CRE-support/sample.bed | awk '$8 == "transcript" {print $0}' > ../NO-CRE-support/RefSeq-knownGene-transcriptonly.bed    

less -S ../NO-CRE-support/non-coding-NO-CRE-superfile.tab | awk '{print $1}' |
fgrep -w -f - ../true-unknown-noncoding.sorted.bed | less -S | awk '$8 ==
"transcript" {print $0}' > ../NO-CRE-support/non-coding-NO-CRE-transcripts.bed     

awk '$3 == "transcript" {print $0}' $REFBED | bedtools intersect -a non-coding-NO-CRE-transcripts.bed -b stdin -wb   

awk '$3 == "transcript" {print $0}' $REFBED | bedtools intersect -a
non-coding-NO-CRE-transcripts.bed -b stdin -wb  | cut -f10,19 | awk '{print $2,
$12}'| sort -u -k1,1 | sed 's/["";]//g' >
./RefSeq-Overlap/transcript-Ids-gene-Ids-overlap.tab
```
* 45 transcripts which are annotated in the RefSeq RefGene catalogue. Remove these from the non-coding no-cre set.     

`cut -d" " -f1 transcript-Ids-gene-Ids-overlap.tab | fgrep -v -w -f - ../non-coding-NO-CRE-superfile.tab > ../non-coding-NO-CRE-superfile2.tab`
* Renamed it and removed the "2" from the end of the name. 
* 1048 transcripts to play with now.     

**Need to perform the identical process with the protein coding set.**     


#### Separate into Distal Enhancer-like Transcripts 
* These transcripts contain a ENCODE-cis-regulatory mark associated with distal enhancer activity. 
* I will be able to profile them further once they are divisioned. For example,
the presence of short conserved DNA regions, or the enrichment of
particular transcriptional factor binding sites (neuronal associated?) 
* Look at their coverage levels, length, exonic content. **Intersection with in house datasets will be valuable**
* GRO-seq, net-seq datasets intersection to conclusively say these are enhancer like transcripts rather than bona-fide lncRNAs.     

`awk '$8 == "transcript" {print $0}' $NCBED | bedtools intersect -a stdin -b local.bed -wb | less -S | cut -f
10,20 | cut -d ";" -f 1,6 | awk '$3 ~ '/^dELS/' {print $0}' | uniq | awk '{print $2}' | sed 's/["";]//g' > non-coding-dELS-intersect.lst`     
- 419 transcripts, - some have additional CTCF binding sites along with the enhancer signatures.     

Extract the corresponding entries from the Superfile. 

`fgrep -w -f non-coding-dELS-intersect.lst ../tmp-newest-superfile.lst > non-coding-dELS-superfile.tab`      

#### Separate into Proximal Enhancer-like Transcripts and Promoter like Transcripts

* I have noticed that many of the promoter-like signatures (PLS) also contain a
proximal enhancer like signature (pELS) - as well as varying amounts of other
overlaps such as CTCF, DNASE H3H4me3 binding. 
* To parse out an enhancer from a PROMPT, the directionality, and distance to
the neighbourng transcript will be used. If a transcript is within 5kb (max),
and is being transcribed in the opposite direction to the adjacent gene, it is
likely a PROMPT. Further support for this if it does not have CTCF binding. A proximal enhancer should have the same directionality.   

```
awk '$8 == "transcript" {print $0}' $NCBED | bedtools intersect -a stdin -b
local.bed -wb | less -S | cut -f 10,20 | cut -d ";" -f 1,6 | awk '$3 != "dELS"
{print $0}' | awk '$3 != "dELS,CTCF-bound" {print $0}' | awk '$3 !=
"CTCF-only,CTCF-bound" {print $0}' | awk '$3 != "DNase-H3K4me3,CTCF-bound"
{print $0}' | sort | uniq | awk '$3 != "DNase-H3K4me3" {print $0}' >
non-coding-pELS-PLS-intersect.lst   
```    
* Removed overlapping hits from dELS dataset - so only unique hits here. 
* 328 hits from 196 transcripts - some have overlapping signatures and will be addressed further on.    


#### CTCF bound, DNASE H3K3me3 accessible transcripts 
* These transcripts here will consist of both promoter-associated, and
enhancer-associated transcripts. They will be parsed and divided into the
respective subsets based upon their directionality, and distance to the
neighbouring gene.    

DNASE-H3K4me3 binding only 
```
awk '$8 == "transcript" {print $0}' $NCBED | bedtools intersect -a stdin -b
local.bed -wb | less -S | cut -f 10,20 | cut -d ";" -f 1,6 | awk '$3 != "dELS"
{print $0}' | awk '$3 != "dELS,CTCF-bound" {print $0}' | awk '$3 !=
"CTCF-only,CTCF-bound" {print $0}' | awk '$3 != "DNase-H3K4me3,CTCF-bound"
{print $0}' | sort | uniq | sort -u -k2,2 | awk '$3 == "DNase-H3K4me3" {print
$0}' > non-coding-DNaseH3K4me3-CTCF-intersect.lst   

# CTCF only-CTCF bound
awk '$8 == "transcript" {print $0}' $NCBED | bedtools intersect -a stdin -b
local.bed -wb | less -S | cut -f 10,20 | cut -d ";" -f 1,6 | awk '$3 != "dELS"
{print $0}' | awk '$3 != "dELS,CTCF-bound" {print $0}' | sort -u -k2,2 | awk
'$3 = "CTCF-only,CTCF-bound" {print $0}' >>
non-coding-DNaseH3K4me3-CTCF-intersect.lst
```  
* Remove those overlapping with distal enhancer set.    
* 280 hits - now to break these down into the PLS/pELS, DELS sets.    


* **Now that the transcripts which contain CRE signatures have been underdone a
  first pass filtration, the transcripts within each set will be further
  scrutinated and separted into their appropriate sets based on 1)
  directionality, and 2) distance to closest gene.**     

Attempt to assign these DNASE and CTCF bound transcripts to the PROMPTS and
enhancer categories based on their sequence features, such as directionality
and distance to closest gene.      

How many are transcribed in the opposite direction to the closest gene?   
`less -S non-coding-DNaseH3K4me3-CTCF-superfile.tab | awk '$5 != $8 {print $0}' | wc -l`    
* 233 - most of them.      

Of these, how many are less than or equal to 3kb away from the closest gene? Using 3kb to stay consistent with definitions.   
`less -S non-coding-DNaseH3K4me3-CTCF-superfile.tab | awk '$5 != $8 {print $0}'
|  awk  'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($9) < 3000)
print $0}' >
./H3K4me3-DNAse-CTCF-bound/H3K4me3-CTCF-bound-Potential-PROMPTS.tab`    
* 216 - the majority. Increasing the threshold to 4kb only included 1 additional transcript, so the threshold was kept at 3kb.    
* The remainder are perhaps distantly transcribed PROMPTS (less likely),
enhancers, pseudogenes, or hopefully, lincRNAs with independent promoters. They
have been output in a separate file as seen below.     

>  …H3K4me3-CTCF-bound-Divergent-to-closest.tab 

Those with a distal enhancer like signature but display PROMPT positional features.    
`less non-coding-dELS-intersect.lst | fgrep -w -f -
../H3K4me3-DNAse-CTCF-bound/H3K4me3-CTCF-bound-Potential-PROMPTS.tab >
../H3lK4me3-DNAse-CTCF-bound/H3K4me3-CTCF-bound-Potential-PROMPTS-withdELS.tab`        



From this set, extract those which have the sequence characteristics of proximal, and distal enhancers.

Proximal enhancers are those which are transcribed from the same strand as the closest gene, and are located less than or equal to 3kb upstream.    

```    
less -S non-coding-DNaseH3K4me3-CTCF-superfile.tab | awk '$5 == $8 {print $0}'
| awk '$5 == "+" {print $0}' | awk '$9 > 0 {print $0}' | awk '$9 < 3000 {print
$0}' >
./H3K4me3-DNAse-CTCF-bound/non-coding-DNaseH3K4me3-CTCF-proximal-enhancer-like.tab

less -S non-coding-DNaseH3K4me3-CTCF-superfile.tab | awk '$5 == $8 {print $0}'
| awk '$5 == "-" {print $0}' | awk '$9 < 0 {print $0}' | awk '$9 > -3000 {p
rint $0}' >>
H3K4me3-DNAse-CTCF-bound/non-coding-DNaseH3K4me3-CTCF-proximal-enhancer-like.tab
```
* 17 hits here 

Distal enhancers are those which have the same criterea as proximal enhancers,
but are located further than 3kb upstream the closest gene. DeSanta uses 20kb
away as a criteria for an enhancer, so be weary of transripts located VERY upstream.      

```     
less -S non-coding-DNaseH3K4me3-CTCF-superfile.tab | awk '$5 == $8
{print $0}' | awk '$5 == "+" {print $0}' | awk '$9 > 0 {print $0}' | awk '$9 <
3000 {print $0}' >
./H3K4me3-DNAse-CTCF-bound/non-coding-DNaseH3K4me3-CTCF-distal-enhancer-like.tab

less -S non-coding-DNaseH3K4me3-CTCF-superfile.tab | awk '$5 == $8 {print $0}'
| awk '$5 == "-" {print $0}' | awk '$9 < 0 {print $0}' | awk '$9 < -3000 {print
$0}' >>
H3K4me3-DNAse-CTCF-bound/non-coding-DNaseH3K4me3-CTCF-distal-enhancer-like.tab
```
* Removed overlap with proximal enhancer like
* 3 hits here 


Create a sub-set which may be Polymerase run-on/rippling of perhaps 3' UTR extensions. Escaping XRN2 cleavage.          

### Parse out likely proximal enhancers and PROMPTS 
* From the PLS-pELS subset, the transcripts which are transcribed in the same
direction as the closest gene, and are located upstream to them, these are most
likely to be proximal enhancers, or even potentially transcripts coming off of
promoters themselves (unlikely).   

```   
awk '{print $2}' non-coding-pELS-PLS-intersect.lst |sort | uniq |  sed
's/["";]//g' | fgrep -w -f - ../tmp-newest-superfile.lst | awk '$5 == $8 {print
$0}' | awk '$9 > 0 {print $0}' > non-coding-pELS-highconf.tab     

awk '{print $2}' non-coding-pELS-PLS-intersect.lst |sort | uniq |  sed
's/["";]//g' | fgrep -w -f - ../tmp-newest-superfile.lst | awk '$5 == $8 {print
$0}' | awk '$9 < 0 {print $0}' >> non-coding-pELS-highconf.tab    
```  
* Remove those overlapping dELS already.     
* 14 unique transcripts.    
* This suggests that the remainder are most likely PROMPTS.    

Divide the PROMPTS-PLS into their own subset.   
`awk '{print $1}' non-coding-pELS-highconf.tab | fgrep -w -v -f -
../non-coding-pELS-PLS-superfile.tab >
../PROMPTS-PLS/non-coding-PROMPTS-PLS-highconf-superfile.tab`    
* 202 transcripts. 
* Vast majority are monoexonic = 172 - 17 = 2, 9 = 3, 3 = 4, 1 = 5. Why these ones are spliced correctly is an interesting question.    
* The furthest distance away from the neighbour gene is 6.5kb, which appears as
an outlier, as the next furthest is 1.6kb, steadily trending downwards in size
at a rather constant rate of 10-30bp.     
* Further weight could be added to this evidence by overlaying with CpG islands? May not be necessary though.     



From the transcripts which do not have a CRE signature of any kind, those with
close similarity to the proximal enhancer signature will be manually inspected
and added to this set if they fit the criteria - upstream of neighbouring gene,
within 2kb (further proximal enhancer in this set is 1788).    

```    
less -S non-coding-NO-CRE-superfile.tab | awk '$5 == $8 {print $0}' | awk
'$5 == "+" {print $0}' | awk '$9 > 0 {print $0}' | awk '$9 <= 2000 {print $0}'
> Potential-Prox-Enhancers.tab     

less -S ../non-coding-NO-CRE-superfile.tab | awk '$5 == $8 {print $0}' | awk
'$5 == "-" {print $0}' | awk '$9 < 0 {print $0}' | awk '$9 >= -2000 {print $0}'
>> Potential-Prox-Enhancers.tab     
```
* 36 potential proximal enhancers. Likely possible to parse manually.     
* The majority are under 3kb from the gene downstream. Similar distribution
compared to the bona-fide proximal enhancer transcripts. 
* Compare to bonafide-enhancer e.g. exonic number, length.    


Seperate out a potential distal enhancer set. This will be more challenging as
they may comprise legitimate lincRNAs and other ncRNA varieties, aswell as
pseudogenes. This can be addressed further later on, so as a starting point,
everything between 20-2kb upstream a gene on the same strand will be classed as
a potential distal enhancer.     

```    
less -S non-coding-NO-CRE-superfile.tab | awk '$5 == $8 {print $0}' | awk '$5
== "-" {print $0}' | awk '$9 < 0 {print $0}' | awk '$9 >= -20000 {print $0}' |
awk '$9 < -2000 {print $0}' >>
./Potential-Enhancers/Potential-Distal-Enhancers-Superfile.tab

less -S non-coding-NO-CRE-superfile.tab | awk '$5 == $8 {print $0}' | awk '$5
== "+" {print $0}' | awk '$9 > 0 {print $0}' | awk '$9 <= 20000 {print $0}' |
awk '$9 > 2000 {print $0}' >>
./Potential-Enhancers/Potential-Distal-Enhancers-Superfile.tab
```     
* 93 potential distance enhancers of 20kb away

#### Parse out likely PROMPTS+Promoter associated transcripts which lack a CRE signature         
* It is my belief that these transcripts can best be differentiated based upon
their distance to the closest gene, and their directionlity. PROMPTS arise very
closesly (most under 1kb it seems), in an oppposite direction of transcription
to the closest neighbour gene, and possibly on protein-coding genes.      
* I will use this approach to add the transcripts which mirror this pattern,
but don't have an existing CRE signature associated with them.      

Sample code for inspection     
` less -S non-coding-pELS-PLS-superfile.tab | awk '$5 != $8 {print $0}' | awk
'{print $1,$9,$6}' | less -S | sort -k2,2 -r -Vd | less -S   `


Remove the potential proximal enhancers that have been extracted from the
original combined set. The potential PLS-PROMPTS will be output into their own
directory, containing a .bed file of *potential* transcripts which may overlap
CpG signatures.    

`awk '{print $1'} Potential-Prox-Enhancers.tab | fgrep -w -v -f -
../non-coding-NO-CRE-transcripts.bed >
../PLS-PROMPTS/Potential-PROMPTS-PLS-transcripts.bed`     


Intersect with CpG islands, those transcripts without any CRE signature. 

```   
CPG=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/Cpg-islands/mm10-cpg-islands-sorted.bed   

bedtools intersect -a Potential-PROMPTS-PLS-transcripts.bed -b $CPG -wa | head     
```    
* No overlap at all, meaning there is minimal support for these via CpG islands.
* Now, look at distance and directionality.    

Next Stage.
Extract transcripts which have the opposite directionality to the closest gene.    

`less -S Potential-Prox-Enhancers-Superfile.tab | awk '$5 != $8 {print $0}' - | wc -l`     
* 421 candidates.    

As per several PROMPT papers, I will be capping the designation of PROMPTS as
3kb divergent and upstream - even though they are typically seen withint 200nt
(nucleosome window it seems). This assessment has also been made with the input
of the data from the high confidence PROMPTS in the other sets (those with CRE
signatures). This is a difficult decision to make, as increasing this border to
4kb includes more transcripts, and so on upwards. Perhaps as a counterbalanced
strategy, those transcripts intiating further than 3kb, but less than say x kb,
can be analyzed to look for enrichment of PAS, and U1 site depletion, another
feature characteristic of PROMPTS, which may also allow them to be lumped into
this set? This compromise may be wise, but may also dip into the "intergenic"
threshold. What about the length of their transcripts? This can be taken into
consideration, as some PROMPT papers comment that some very long transcripts
are produced upon exosome/quality control depletion... Is there a sweet spot in
the distance from the closet gene, and the potential for it being designated a
prompts? e.g. after 5kb do they then begin losing the PAS enrichment etc.? Are
there particular discirimators to aid in the annotation process.       

`less -S Potential-Prox-Enhancers-Superfile.tab | awk '$5 != $8 {print $0}' - |
awk 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($9) < 3000) print
$0}' > ../PLS-PROMPTS/high-conf-potential-PROMPTS-PLS-3kb-Superfile.tab`    
* 120 transcripts - now what to do with the remaining ~300, some of which may also resemble a similar PROMPT like signature.     
* The overwhelming majority are monoexonic! Booya. 108 = 1, 9 = 2, 3 = 3.      
* The ones which are spliced, why are they escaped?     

Output prompts further than 3kb, but less than 5kb upstream.  
`less -S ../Potential-Enhancers/Potential-Prox-Enhancers-Superfile.tab | awk
'$5 != $8 {print $0}' - | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {i
f (abs($9) > 3000) print $0}' | awk 'function abs(x){return ((x < 0.0) ? -x :
x)} {if (abs($9) <= 5000) print $0}' >
low-conf-distnt-PROMPTS-PLS-Superfile.tab`     
* 58 transcripts emerging between from 3-5kb upstream. Some of these may in fact be lncRNAs or other elements.      


#### So what are the remainder of the transcripts that don't resemble enhancers and PROMPTS?     

Created a set of transcripts which are NOT part of any of the above described sets of genetic elements.
Within this, there will be pseudogenes, 3' Pol-II run-on/Extensions, and perhaps even bonafide non-coding genes.    

```    
awk '{print $1}' low-conf-distnt-PROMPTS-PLS-Superfile.tab >  ../not-lncRNAs.lst   
awk '{print $1}' high-conf-potential-PROMPTS-PLS-3kb-Superfile.tab >>  ../not-lncRNAs.lst   
awk '{print $1}' Potential-Prox-Enhancers.tab >>  ../not-lncRNAs.lst   
awk '{print $1}' Potential-Distal-Enhancers-Superfile.tab >>  ../not-lncRNAs.lst     
``` 
* 307 transcripts belong to the above mentioned sets - the remainder are free lunch.     


#### Subset of 3' pol-II run on, 3' UTR extension transcripts.      
* These should be transcribed from within the last x nucleotides of the 3' end
of the closest gene, on the same strand.
* There are 741 potential transcripts to play with.     
* They will be dispursed across the transcripts without a CRE-signature. They
shouldn't be present within the enhancer, or PROMPT datasets, based on both
their directionality, and proximity to the 3' end, rather than the 5' end.     
* DeSanta (2010) used 10kb downstream as a cutoff for polymerase run-on
filtering. However, as I believe the approach that they used was not
RNA-sequencing, but instead Ser5P, Ser2P Chip (Pol-II binding sites), they used
a stricter, more sensitive threshold for Pol-II initiation. 
* I will use anything up to 3kb downstream as pol-II run on and spurious intiation - especially if only very lowly expressed.    

```    
less -S other-genetic-elements-superfile.tab | awk '$5 == $8 {print $0}' | awk
'$5 == "+" {print $0}' | awk '$9 < 0 {print $0}' | awk '$9 >= -3000 {print $0}'
> ./3-Prime-Runon-Extensions/3-Prime-runon-extension-Superfile.tab     

less -S other-genetic-elements-superfile.tab | awk '$5 == $8 {print $0}' | awk
'$5 == "-" {print $0}' | awk '$9 > 0 {print $0}' | awk '$9 <= 3000 {print $0}'
>> ./3-Prime-Runon-Extensions/3-Prime-runon-extension-Superfile.tab     
```      
* 79 transcripts are associated with this set. 
* The remainder are thus our lincRNA - ncRNA candidates, and pseudogenes.    


#### lincRNA - ncRNA set.      
*  These will be any transcript which hasn't been allocated to another of the
other sets described above. As such, they do not resemble enhancer, PROMPTS, 3'
run ons, pseudogenes, and hopefully, represent at least some bona-fide lncRNA genes.    

`awk '{print $1}' 3-Prime-runon-extension-Superfile.tab | fgrep -w -v -f -
../other-genetic-elements-superfile.tab >
../lincRNA-ncRNA/lincRNA-non-coding-set-superfile.tab`   
* 662 potential lincRNAs - keep in mind, I have not filtered these for pseudogenes yet.     
* Are there any differences in this set compared to the others? Exon count? Length? Expression/Coverage?     

#### BLASTN on non-coding set - Pseudogene search.      
Run blastn against the nucleotide sequences of the non-coding set as it may be
the surest way to detect high hits to pseudogenes.    
```
REFSEQ=/media/labpc/Disk-2/Bioinformatics-Computational/blastdb/refseq/knownGene/knownGene        

blastn -db $REFSEQ -query ../*.fa -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > true-unknown-noncoding_blastn.outfmt
``` 
* 778 positive matches.     
* List coverage levels - many of these are barely detectable it appears - the
combination of very low expression, paired with high sequence similarity, has
contributed to their unknown status. Does this have implications for
expression mismeasurement? Are the mappers falsely assigning reads to their
paralogue?     
* Has turned into a gold-mine for pseudogene discovery - verifying my initial
beliefs in that performing a blastn and parsing for highest alignment lengths
may be a viable tool - Almost all have been pseudogenes.    
* Have to remove positive ref-seq intersects from this set.    

### Split both non-coding and protein coding candidates based on coverage level
* High = 100 and over, Mid = 25 - 100, Low = < 25.    

Protein coding set.    
```   
less -S Superfile_Cov_finale.tab | awk '$3 > 100.00 {print $0}' - > Superfile_Cov_high_finale.tab   

less -S Superfile_Cov_finale.tab | awk '$3 < 100.00 {print $0}' - | awk '$3 > 25.00 {print $0}' - > Superfile_Cov_mid_finale.tab    

less -S Superfile_Cov_finale.tab | awk '$3 < 100.00 {print $0}' - | awk '$3 < 25.00 {print $0}' - > Superfile_Cov_low_finale.tab      
```  
* 5 transcripts (some with multiple ORF predictions) in the high set.    
* 17 transcripts in the mid set.   
* 769 transcripts in the low set.     
* GOOD LUCK PARSING THESE SONNY!    

Non-coding set.   
```   
less -S tmp-new-superfile.lst | awk '$2 > 100.00 {print $0}' - > NC_Superfile_Cov_high.lst

less -S tmp-new-superfile.lst | awk '$2 < 100.00 {print $0}' - | awk '$2 > 25.00 {print $0}' - > NC_Superfile_Cov_mid.lst    

less -S tmp-new-superfile.lst | awk '$2 < 100.00 {print $0}' - | awk '$2 < 25.00 {print $0}' - > NC_Superfile_Cov_low.lst     
```     
* 7 transcripts in the high set.    
* 30 in the mid set.    
* 1755 in the low set.    

### Small RNA scans. 
* I am most interested in scanning for small RNAs from within the non-coding candidates (both with and without CRE signatures).     
* I will also scan the potential protein coding candidates also.    
* snoRNAs (HA Box, and CD) using SNOscan and SNOgps (for the different snoRNAs respectively). 
* tRNAs - using tRNA-SE (2.0) by Lowe Lab.   
* Attempt to do SRP scan, but this may be more challenging.     

#### tRNA scan

On the non-coding candidates at first. 
$NCNOS are just all of the non-coding candidates in fasta format.
`tRNAscan-SE -y $NCSNOS -o trna-SE.out -b trna-SE.bed -a trna-SE.fa`
* 

With highest sensitivity
`tRNAscan-SE --max -y $NCSNOS -o trna-SE-hisens.out -b trna-SE-hisens.bed -a trna-SE-hisense.fa`     


#### SRP scan     
* Most of the relevant details can be found on these following links. Much of it is aged, but I think, still relevant.    
* http://bio.lundberg.gu.se/srpscan/rnascan_doc.html
* http://bio.lundberg.gu.se/srpscan/paper.html   
* http://compbio.fmph.uniba.sk/rnarobo/
* https://github.com/rampasek/RNArobo   


`rnabob -F -c Eukaryotes_Others.des $ALLSRP >
/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/non-coding-set/small-RNA-scans/SRP-all-pc-nc-rnabob.out`    


### Pseudogene Parsing. 
* In this section I will endevour to compile a high-confidence pseudogene set
from the various blastn (and potentially tRNA scans) scans.   
* I have already manually identified roughly 20 pseudogenes, but there are many
more potentially hits which need to be addressed.      

#### Remove entries from non-coding blastn search which intersect known RefSeq annotations.      
* Filtration!     

`awk '{print $1}' transcript-Ids-gene-Ids-overlap.tab | fgrep -w -f -v
../../blastn/true-unknown-noncoding_blastn.outfmt >
../../blastn/true-unknown-noncoding_blastn-2.outfmt` 
* 32 hits. Proceed with more clarity now.   
* Copied file from the above directory to non-coding-set/pseudogenes/potential-pseudogenes-blastn.tab


#### Non-coding hit parsing 
* I may have missed some based upon my stringency criteria, but I drew the line
at a particular point based upon my manual searches of low scoring hits -
almost every single manual search of less stringent hit resulted in detection
of a TE. 
* All hits must be more than 95% similar, have a length over 180, and a score
 greater than 200.   
* The set was manually parsed and searched on BLAT.     

`less -S potential-pseudogenes-blastn.tab | sort -k12,12 -rn | awk '$3 > 95.00
{print $0}' | less -S | sort -k3,3 | awk '$4 > 100 {print $0}' | awk '$12 > 200
{print $0}' | awk '$4 > 180 {print $0}' | less -S`    

Manually compiled a tab separated file containing the IDs of the transcript,
with the parent gene in the adjacent column.    
* 61 candidates.       

Create a corresponding .bed file.     
`awk '{print $1}' nc-pseudogenes-hiconf.lst | fgrep -w -f -
../true-unknown-noncoding.sorted.bed | awk '$8 == "transcript" {print $0}' >
nc-pseudogenes-hiconf.bed`     


#### Protein-coding candidates hit parsing     
* First I need to remove the transcripts which match a refseq annotation that
hasn't been listed on the gencode annotation that I used.    

```     
less -S transdecoder-candidates-unknown.gtf | awk '$3 == "transcript" {print $0}' > tmp-trans.gtf

awk '$3 == "transcript" {print $0}' $REFGTF | bedtools intersect -a ./tmp-trans.gtf -b stdin -wb > transdecoder-candidates-refseq-intersect.gtf

less -S transdecoder-candidates-refseq-intersect.gtf | cut -f9 | awk '{print $2}' | sed 's/["";]//g'
| sort | uniq > transdecoder-candidates-refseq-intersect-IDs.lst
```      


Remove positive refseq hits from pseudogene candidates.     

```    
reflist=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/TransDecoder/RefSeq-intersect/transdecoder-candidates-refseq-intersect-IDs.lst     

fgrep -v -w -f $reflist true-unknown.fa.transdecoder-blastn.outfmt > true-unknown.fa.transdecoder-blastn2.outfmt
```        

Manually parse hits as done above for non-coding candidates.   
* Created a tab-separated file with each Stringtie ID and corresponding parent gene "pc-candidates-pseudogenes-blastn.lst"    
* 12 pseudogenes.      


#### Intersect pseudogene candidates with short-read RNA seq data (in house).    

Now that both non-coding and protein-coding candidate files have been parsed, we can combine them.   

```    
PCBED=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/TransDecoder/blastn/pseudogenes/pc-candidates-pseudogenes-blastn.bed    


NCPG=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/non-coding-set/pseudogenes/nc-pseudogenes-hiconf.bed
```    

Compute the average coverage/depth of the positive intersects in the short-read
data - answering the question: **if transcription is detected from the
pseudogenic loci in the short-read data, at what sequencing depth has it been detected?**      
* Refer to thread (here)[https://www.biostars.org/p/279140/] for help.       


```
samtools merge $RC $EXT > merged.bam 
# Indexed and sorted 

bedtools bamtobed merged.sorted.bam | sort -k1,1 -k2,2n > out.ssorted.bed 

cat $PCBED $NCPG | sort -k1,1 -k2,2n | bedtools coverage -sorted -a - -b
out.ssorted.bed -mean >
/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Pseudogene-Pred/short-read-support/NC-PC-PS-coverage-mean.bed     

# Calculate average coverage.    
less -S *.bed | awk '{print $NF}' | awk '$1 != 0 {print $0}' | awk '{total+=$1} END {print total/NR}'
```    

* 4.19 coverage at each base for non-coding candidates - so quite low coverage as expected.    
* 11.1659 for combined protein-coding and non-coding candidates  - What is the max, min and medium-mode?    
* Max is 306 (wow? - what?) - min is 0.33. Most are under 6x coverage.        


Intersect novel-pseudogenes with existing merged short-read data.    

* It seems the short-read stringtie.gtf only resulted in 1 intersect, leading me to believe
that it may not be reconstructing the transcripts correctly, and perhaps
assigning them to the parent gene. As such, I will simply convert the mapped
.bam files to .bed format, and intersect them as such, looking for any traces of
"transcriptional activity", and inferring from there onwards.        
* GRO-seq data will also be of help here.    

```     
SRBED=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/pending-uplo
ad/HISAT2/merged/merged/out.ssorted.bed     

cat $PCBED $NCPG | sort -k1,1 -k2,2n | awk '$8 == "transcript" {print $0}' |
bedtools intersect -sorted -a - -b $SRBED -wb >
PC-NC_PS_Short-read-Intersect.bed

# Create a .bam file based on the regions in the dRNA-pseudogene .bed 
less -S *.bed | cut -f1,2,3,4,5,6,7,8,9,10 | samtools view -@ 10 -L -
/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/pending-upload/HISAT2/merged/merged/out.sorted.bam
> PC-NC_PS_Short-read-Intersect.bam      
```       

cat $PCBED $NCPG | sort -k1,1 -k2,2n | awk '$8 == "transcript" {print $0}' |
samtools view -@ 10 -b -L - RC-EXT.bam >
/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Pseudogene-Pred/NC-PC-Pseudogene-set-transcript.bam  


* 44 pseudogenes are supported by the short-read data.      
* Extract the read ID's of the short-read data from the output column and create .fasta sequences from this.     



#### Download UCSC 60-way MultiZ alignments 
* Most streamlined way to do this is with rsync, but with the addition of --include and --exclude to avoid downloading incomplete haplotigs.     
`rsync -amvz --progress
rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/multiz60way/maf/chr*.maf.gz
--exclude='*Un*' --exclude='*random*' ./`

These .maf.gz can be converted to individual .bed files with conservations
scores contained in their own column. This allows one to specify a specific
chromosomal locus and evaluate the conservation score of said region. This is
alignments for the specific transcripts that I am interested in.. ?      


#### Anti-fam scan 
* Search for spurious ORF predictions which has positive PFAM hits - rule out false positives.    

```     
hmmsearch --cut_ga AntiFam.hmm
ANTIFA=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/TransDecoder/blastn/true-unknown-pc-candidates-fullfasta.fa

hmmsearch --domtblout out.out --cut_ga
/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/Anti-Fam/AntiFam.hmm
$ANTIFA    
```  
* No hits - proceed with ORF prediction hits.     



#### Intersect non-coding CRE-distal enhancer candidates with enhancer signatures
* From Greenberg 2010 (Kim) neuronal enhancer dataset. I will intersect with CBP and H3K4me1 in untreated cell and KCL stimulated.  

Sample script below. 

`less -S /Distal-Enhancer-dELS/non-coding-dELS-intersect.lst | fgrep -w -f -
true-unknown-noncoding.sorted.bed | bedtools intersect -a - -b
GSM530204_H3K4Me1_H3K4Me1_KCl_B2_E120.bed -wb | sort -k19,19 -rn | less -S   
`       
* No intersects above baseline/background levels of binding ~1.      
* Perhaps we may have signal with p300, H3K27ac, or even H3K4me3, but for now, my bandwidth has been used up.     



---------------------------------------------------------------------------------------------------------------------------------


Immediate


- UCSC MultiZ Multiple Alignments. Following DAVE tangs tutorial to get the
  tracks up and running. Figure out how to extract multiple alignments for the transcripts I am searching for.             
- Mpe16 Pseudogene/lincRNA RNAcode. - Other potential, high confidence protein
  coding genes into RNAcode?    

- Determine what the next steps are with the pseudogene set: (1) Look for
  intersection in Short-Read data (2) Look for intersection in bambu, (3)
  Attempt reconstruction in simulated data (is it repeat content, similarity to
  existing UTR, low coverage level which has resulted in their poor detection?)
  - (4) Attempt to predict new pseudogenes via 3'UTR + genome which has had
  annotated regions removed.      
- Pseudogene multiple alignments - stop codons, frameshifts etc. Take the
  respective distance upstream of the pseudogene, which corresponds to the
  length of the parent gene, then align it, to see if the remainder of the gene
  is also present - **this can be done to rule out the notion that only the 3'
  UTR is being transposed/duplicated?**            

Protein-coding candidates.      

Final Wrap Up    
- SQANTI for quality control reports and potentially additional features. 


**Potentially bona-fide lncRNAs**   
- Ideally have some enriched conservation regions, splice sites, promoter region etc.     

**Enriched in Low-complexity regions** 
- Download refseq genes.fa (repeat masked lower case) and calculate lower case nucleotide composition of each gene.     
- RepeatMasker regions? Is this enrichment of TE, over or under the 'average' for genes? See Goke benchmark paper.     
- Percentage of the sequence that is repeat masked? A basic calculation is possible surely. Ratio of lower-case to upper-case characters.     
`grep -o [a-z]
/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm24/
| tr -d "\n" | wc -m` 
*WORK IN PROGRESS* 
- The notion being pursued here, would shed light on perhaps one of the reasons
  these transcripts have not been detected with short-read technologies, but
  are now being uncovered with longer reads.       


**Evolutionary Support**   
- Plot PhastCons, PhyloP across gene body - potentially also around splice sites.  
- Evolutionary rate analysis of high confidences genes - under selection? Multiple species alignment needed.    


**Data presentation**   
- Plot Poly-A cluster, then plot poly-a tail from squiggle data.     
- Begin generating simple plots in R to visually represent the above results - this will tie into the presentation/slideshow 


**Pseudogenes**    
- Do they have poly-a tails? - makes sense that only the 3' most side is being
transcribed, afterwhich the transcript is promptly terminated as it lacks the
necessary signals to, or contains negative signals, for the passage of
polymerase 
- Do the 3' most pseudogenes terminate at a similar distance, due to a
degradation signal. Take upstream sequence of the pseudogene, and perform a
hexamer enrichment/scan to identify any potential cleavage/degradation signal?   
- Align the ends of the pseudogenes 
- ONT sequences 3' to 5' and thus, it captures the poly-Adenylated pseudogene fragment before it is potentially degraded? One advantage of ONT.     
