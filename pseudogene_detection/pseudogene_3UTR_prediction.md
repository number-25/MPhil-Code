# Document outline 
Herein, I will attempt to perform an initial pseudogene detection strategy,
based upon the observations within the directRNA dataset (noncoding blastn
paired with manual BLAT). Extracting the most positive BLAST hits, and
searching them manually on BLAT, I uncovered a consistent trend in which the 3'
UTR of ~60  protein coding genes was being detected at other genomic regions
(almost in complete identity). Furthermore, these hits were also frequently
flanked by tranposable elements, and had similarity to mRNA in other species,
indicating some pseudogenic qualities, along with conservation blocks. 

Based upon these trends, I pondered whether it would be possible to predict,
and detect additional unnanoted pseudogenes by extracting the last ~500bp of
every genes 3' UTR, store these in their own .fasta format, extract all the
genomic nucleotides which don't have any genetic annotation ascribed to them,
then map the 3'UTR.fa to this index.      


## Aims 
* Perform prediction described above.
* Discover un-annotated pseudogenes.    



## Analysis 

### Extract the last 500 basepairs of the 3' UTR

Use (seqkit)[https://bioinf.shenwei.me/seqkit/usage/#subseq]. This can also be
done with
(biopython)[https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec19] as demonstrated on (Biostars)[https://www.biostars.org/p/710/].               

`seqkit subseq -r -500:-1 gencode.vM24.transcripts.fa >
/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Pseudogene-Pred/3UTR-scan/gencode.vm24.transcripts-3UTR-last500.fa`     

Merge gencode (vm24) comprehensive, refseq knowngene, and gencode pseudogenic (vm25) and tRNA annotations.     

First, I must assign the annotation to convenient variables. 
```     
GENBED=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm24/filters/gencode.vM24.annotation-tx-only.bed    

GENPS=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm25/gencode.vM25.2wayconspseudos.sorted.bed    
```

Sample code
cat ./local-delete-later/refbed ./local-delete-later/genbed
./local-delete-later/genps | bedtools merge -i - | sort -k1,1 -k2,2n | bedtools
complement -i - -g ./local-delete-later/genome-order | bedtools getfasta -fi
../../mm10/mm10.fa -bed - | less -S | fold -w 60 >
/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Pseudogene-Pred/3UTR-scan/extra-genic-genome.fa

As the strategy of merging the annotation files is not proving fruitful, and is
leading to errors which I am unable to resolve, I will simply use a single
annotation (gencode.vm24) to create the extragenic.fasta. After mapping the
3'UTR regions to this extra-genic region, I will convert from .bam to .bed
file, and perform an intersect (-v) against the RefSeq annotation, and the
pseudogenic annotation.      


Complement the gencode vm24 transcriptome annotation (comprehensive), then take
the nucleotide sequences of these extra-genic regions and output them into
their own .fasta file.   
* Had to sort the genome size file with k1,1 > genome-order so that the analysis would proceed correctly.     


```   
bedtools complement -i ./local-delete-later/genbed -g
./local-delete-later/genome-order | bedtools getfasta -fi ../../mm10/mm10.fa
-bed - | fold -w 60 >
/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Pseudogene-Pred/3UTR-scan/extra-genic-genome.fa
```   

Remove sequences shorter than 30nt with Biopython (script)[https://biopython.org/wiki/Sequence_Cleaner]    
`python Sequence-cleaner.py gencode.vm24.transcripts-3UTR-last500.fa 30 30`      

### Map with minimap2 

Map with minimap2 (version 2.22)     
`minimap2 -t 8 -a 3UTR-last500-index.mmi extra-genic-genome.fa > trial.sam`
* According to samtools flagstat, there are no successfull alignments at all.... hmmmmm


### Install local BLAT

* I followed this guide to install BLAT correctly - https://bioinformaticsonline.com/pages/view/37677/installing-blat-on-linux
* https://genome-test.gi.ucsc.edu/~kent/src/  

```     
echo $MACHTYPE   
MACHTYPE=x86_64   
export MACHTYPE   
echo $MACHTYPE   
vim ~/.bashrc    # export PATH=~/bin/x86_64::$PATH   
source ~/.bashrc   
mkdir ~/bin/$MACHTYPE   
vim inc/common.mk   # BINDIR=/usr/local/bin   
sudo make   
``` 

Run local BLAT with options. Useful reading
(here)[https://genome.soe.ucsc.narkive.com/vgbggUcP/blat-11-ooc-file-from-blatsuite-zip]         
`blat -t=dna -q=dna -dots=1000 extra-genic-genome.fa clear_gencode.vm24.transcripts-3UTR-last500.fa output.psl`      
* After perusing through the output, I observed many, many matches for each
transcript - low, high scoring, the whole lot. As I wanted only the top scoring
hits, I changed the "minScore" and "minIdentity" parameters, and reran the
analysis.    

Minimum score of 400, and a minIdentity of 95%.    
`blat -t=dna -q=dna -dots=1000 -minScore=400 -minIdentity=95 extra-genic-genome.fa clear_gencode.vm24.transcripts-3UTR-last500.fa output.psl`    
* Needs more refining - far too many hits.  


Mask lower case (repeats), use recommended "fastMap" for DNA to DNA searches, and use newly created .2bit database for speed. 

`blat -fastMap -mask=lower -t=dna -q=dna -dots=1000 -minScore=400
-minIdentity=95 extra-genic-genome.2bit
clear_gencode.vm24.transcripts-3UTR-last500.fa output.psl`   
* Much quicker, more specific with some hopeful candidates.     

Remove entries mapping to incomplete and random chr contigs, more than 50 repeat masked matches, and inserts greater than 20bp.    
`less -S output-noheader.psl | grep -v "random" | grep -v "Un" | awk '$3 < 50 {print $0}' | awk '$8 < 20 {print $0}' > output-noheader2.psl`   

See if I can convert this to .bed format for easier parsing and future intersection.   
https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/psl2bed.html    

`cat output.psl | psl2bed > output.bed`    

Use kent-utils program pslToBed 

`~/Bioinformatics-Programs/UCSC-utils/pslToBed output2.psl output2.bed`      
* Doesn't seem to be formatted correctly - columns are spaced strangely.      
* May have to manually create .bed file if the intersections do not work.      


Requires some more manual formatting as the first three columns are not
correct. Columns 4 and 5 are also unimportant.      
`less -S output.bed | sed
's/:/ /g' - | awk '{gsub(/[-]/, "\t", $2)} 1' OFS='\t' - | awk '{print
$1,$2,$3,$6,$7,$8,$9,$10,$11,$12,$13,$14}' OFS="\t"`     
* This brings us to 6995 unique positive hits across the extra-genic genome.    
* Now I will intersect the newly formatted .bed file with additional genome annotations.     


### Intersect with existing annotations.    
* In order to evaluate whether this mapping strategy is relatively effective or
not, I will create a .bed file of the BLAT output, and intersect it with the
gencode pseudogene annotation to determine if it is detecting existing
pseudogenes.    
* Next I will extend the reads by 1250bp upstream and perform another
intersect, to determine whether it can be used to create more complete
pseudogene models which represent homology to known bona-fide genes.     
* Everything outside of these intersections is FREE lunch.      



Intersect against GENCODE vm25 2way pseudogene annotation.     
```   
PCBED=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm25/gencode.vM25.2wayconspseudos.sorted.bed    

bedtools intersect -a BLAT-output-sorted-formatted.bed -b
/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm25/gencode.vM25.2wayconspseudos.sorted.bed
-wb | awk '{print $13,$14,$15,$41} | sort | uniq | wc -l    
```     
* 526 unique overlaps to existing pseudogene annotations - of the 6995 unique BLAT hits.    
* Next, I will take the positive hits which didn't overlap an existing
annotation, extending them by 1200bp in each direction, then perform this intersect
again, to see if I am detecting the 3' head of existing pseudogenes, providing
more complete pseudogene annotations if this fragement so happened to be
transcribed and detected in a sequencing study.       


```    
bedtools slop -b 1200 -i
./Pseudogene-Gencode-Intersect/PC+NC-NO-intersect.bed -g $mm10 | bedtools
intersect -a - -b $PCBED -wb | awk '{print $13,$14,$15,$41}' >
Pseudogene-Gencode-Intersect/PC+NC-pseudogene-slop-gencode-intersect-uniq.tab      

bedtools slop -b 1200 -i
./Pseudogene-Gencode-Intersect/PC+NC-NO-intersect.bed -g $mm10 | bedtools
intersect -a - -b $PCBED -wb >
Pseudogene-Gencode-Intersect/PC+NC-pseudogene-slop-gencode-intersect.bed       

./Pseudogene-Gencode-Intersect/PC+NC-NO-intersect.bed -g $mm10 | bedtools
intersect -a - -b $PCBED -wb | awk '{print $13,$14,$15,$41}' | sort | uniq | wc
-l     

bedtools slop -b 1200 -i PC+NC-NO-intersect.bed -g $mm10 | bedtools intersect -v -a - -b $PCBED > PC+NC-NO-SLOP-intersect.bed

```    
* 209 additional unique hits to known 2way pseudogene annotations. There is a very,
very marginal increase in pseudogene intersection by increasing the slop from
200 to 2000bp - as such, I settled for 1200bp.       


Now intersect with the UCSC-Gencode main pseudogene annotation file. 

```   
UCSCPS=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/GENCODE/vm25/vm25-gencodepseudogenes-sorted.bed    

bedtools intersect -a PC+NC-NO-SLOP-intersect.bed -b $UCSCPS -wb > PC+NC-pseudogene-slop-gencode-ALL-intersect.bed 
```
* 342 additional unique hits to known GENCODE UCSC pseudogenes.     

How many of the remaining hits do we have in the data which do not intersect the pseudogene annotation? 
* 1062 
* Use these as the high-confidence set once ref-seq intersections are performed to remove existing annotations.    
* Remove the 1200bp that we slopped on earlier. 


Perform the final intersect to ensure that these are not mapping to known RefSeq genes.    
```    
bedtools intersect -v -a PC+NC-NO-SLOP-intersect.bed -b $UCSCPS | bedtools slop
-b -1200 -i - -g $mm10 | bedtools intersect -v -a - -b
/media/labpc/Disk-2/Bioinformatics-Co
mputational/Reference-Genomes/RefSeq/mm10/mm10.refGene-transcripts.bed >
../high-conf-unique/PC+NC-no-annotations-hiconf-final.bed      
```    
* Left with 1781 unique BLAT hits which do not intersect Pseudogene, nor RefSeq annotations.     


Many of the manual BLAT searches are still showing hits just downstream of
existing pseudogenes, so I am thinking a double slop strategy is most
appropriate to eliminate this. 
- **UPDATE** turns out that I was not using the complete pseudogene annotation from UCSC, only the 2way predicted pseudogenes.     


Intersect with Repeat Masker files to reduce amount of transcripts mapping to TEs. 

`less -S PC+NC-no-annotations-hiconf-final.bed | bedtools intersect -v -f 0.90
-a - -b
/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/RepeatMasker/mm10/UCSC-RepeatMasker-locations_mm10.bed
| sort | uniq | cut -d"|" -f5 | sort | uniq | wc -l`     
* 137 potential transcripts mapping to unique, extragenic (until now) regions.    

```     
less -S PC+NC-no-annotations-hiconf-final.bed | bedtools intersect -v -f 0.90
-a - -b
/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/RepeatMasker/mm10/UCSC-RepeatMasker-locations_mm10.bed
| sort | uniq  > PC+NC-no-annotations-repeatmasker-hiconf-final.bed    

less -S PC+NC-no-annotations-hiconf-final.bed | bedtools intersect -v -f 0.90
-a - -b
/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/RepeatMasker/mm10/UCSC-RepeatMasker-locations_mm10.bed
| sort | uniq | cut -d"|" -f5 | sort | uniq >
PC+NC-no-annotations-repeatmasker-hiconf-final.lst     
```       


Manually parse the list of transcripts and check with online BLAT.     
* Still see a majority of these mapping to transposable elements despite the repeat masker overlap - perhaps I need to tweak the bedtools intersect.     

Potentially some newly discovered Pseudogenes/duplications? - All of these show very similar
patterns to those discovered in our dRNA with overlap of conservation LOD
blocks.     
* Kif5b-201 on chr2   +   109059034 109059515
* 1110002L01Rik on chr14  +    24392049  24392548
* Foxk2-201 on chr4   +    23230054  23230542  and Chr9, Chr12
* Zfp106-201 on chr14  +    80312603  80313089 and  chr2   +   164306962 164307454  
* Trim33-204 on chr3   -    36116466  36116954     
* Snx4-203 on chr17  -    22094887  22095377      
* Slc5a8-201 on chr10  +    88935561  88936060
* Rgs3-202 on chr19  -    53113808  53114274
* Ptpn9-201 on chr15  +    14621558  14622023    
* Nemf-209 on chr2   +    92811227  92811683 - overlaps a CRE enhancer locus 
* Mthfs-204 on chr9   +    88965122  88965621

Format these into a .bed file manually.

`vim PC+NC-predicted-PS-manual-complete.bed`     

## What's next?
- Overlap these with conservation blocks: this will complement original dRNA
pseudogene conservation block overlap to show that they have similar qualities.   
- Perhaps also see if these loci are transcribed in the GRO-seq data.     
- May need to tweak the BLAT search to avoid matching solely to repeat-mased
regions - the task would be to effectively mask the gencode-transcripts.fa as
it is all upper-case.     

