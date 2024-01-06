# Document outline 
Herein, I will intersect the peptide coding candidates predicted by
TransDecoder with the 60way PhyloP conservation data, in order to look for
transcripts with elevated conservation signal which may be indicative of
protein function.    

This script adds a single exon entry spanning the entire transcripts below all
monoexonic transcript .gtf entries. This allows the phyloP intersection to be
performed on transcrips that may be bonafide peptide coding yet monoexonic, and
restrict the conservation search to the actual exons of potential genes which
are presumably encoding some protein domains.     

```
less -S Superfile_Cov_finale-3.tab | cut -d" " -f5 | sort | uniq | fgrep -w -f
-
../../../StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/true-unknow
n.gtf | ~/Bioinformatics-Programs/gffread/gffread --gene2exon - -T | awk '$3 ==
"exon" {print $0}' | gtf2bed - > true-unknown.fa.transdecoder-exons_FULL-3.bed    
i

less -S true-unknown.fa.transdecoder-exons_FULL-3.bed | bedtools intersect
-sorted -a - -b $PHYLOP -wb >
./PhyloP-Conservation/PC-candidates-transdecoder-PhyloP.bed     
```

`awk '{print $11, $17}' PC-candidates-transdecoder-PhyloP.bed | awk '$2 >= 3 {print $0}' | cut -d" " -f 1 | sort | uniq | wc -l`
* 236 transcripts within this category.      

A great challenge here is validating whether an elevated conservation signal is
distributed over the predicted ORF/peptide of transdecoder, and not within
say, the introns of the gene or outside the potential frame. As such, inspecing the transdecoder gff or bed outputs may prove valuable. 

* It appears that many of the highly-conserved transcripts are overlapping CRE
Enhancer signatures - so it may be informative to intersect these with the
data to narrow down how many are enhancer related.     
* Because of this, may be interesting also to eventually do a meme/hexamer enrichment analysis.     


In a slight detour, chasing understanding and clarity, I have decided that
intersecting the peptide-coding candidates with the ENCODE regulatory data may
be the most suitable. Those transcripts with high conservation and positive CRE
regulatory marks are most likely not functioning as primary transcripts, but
instead, regulatory hubs/platforms for enhancers and so forth. The transcripts
OUTSIDE of this set may be most likely to represent bona-fide peptide coding
elements.    

```
CRE=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Genomes/ENCODE-CRE/encodeCcreCombined.bed 

less -S Superfile_Cov_finale-3.tab | cut -d" " -f5 | sort | uniq | fgrep -w -f
-
../../../StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/true-unknown.gtf
| awk '$3 == "transcript" {print $0}' | gtf2bed - | bedtools intersect -a - -b
$CRE -wb > CRE-intersect/PC-candidates-transdecoder-CRE-support.bed

less -S Superfile_Cov_finale-3.tab | cut -d" " -f5 | sort | uniq | fgrep -w -f
-
../../../StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/true-unknown.gtf
| awk '$3 == "transcript" {print $0}' | gtf2bed - | bedtools intersect -a - -b
$CRE -wb | less -S | awk '{print $11}' | sort | uniq | wc -l    
```     
* 379 transcripts have positive intersects.. wow! 
* So what fraction of these have conservation peaks above 3? Let's find out. 

```    
less -S Superfile_Cov_finale-3.tab | cut -d" " -f5 | sort | uniq | fgrep -w -f
-
../../../StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/true-unknown.gtf
| awk '$3 == "transcript" {print $0}' | gtf2bed - | bedtools intersect -a - -b
$CRE -wb | less -S | awk '{print $11}' | sort | uniq | sed 's/["";]//g' | fgrep
-w -f -
PhyloP-Conservation/PC-candidates-transdecoder-PhyloP-Over3-Superfile.tab >
PhyloP-Conservation/PC-candidates-transdecoder-PhyloP-Over3-CRE-Superfile.tab
```

* 194 do, 2/3's of the set. 115 do not have an intersect.   
* `awk '{print $11, $17}' | sort -k2,2 -rn |` to parse transcript ID and conservation score.     

```less -S Superfile_Cov_finale-3.tab | cut -d" " -f5 | sort | uniq | fgrep -w -f
-
../../../StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/true-unknown.gtf
| awk '$3 == "transcript" {print $0}' | gtf2bed - | bedtools intersect -a - -b
$CRE -wb | less -S | awk '{print $11}' | sort | uniq | sed 's/["";]//g' | fgrep
-w -v -f - PhyloP-Conservation/PC-candidates-transdecoder-PhyloP.bed >
PhyloP-Conservation/PC-candidates-transdecoder-PhyloP-Over3-NOCRE-Superfile.tab```      


* A few examples
STRG.30184.1 - Typical PROMPT signature, with very high conservation, multiple enhancer overlap, mRNA overlap.    
STRG.21636.1 - Transcription through enhancers 
STRG.10101.1 - Another example of transcription through enhancers 
STRG.5217.1  - As above.
STRG.3518.1 - A structural element?     
STRG.15955.1 - A pseudogene? Detected already?   
STRG.19052.1 - A pseudogene? 60x coverage, pfam, uniref hits, 
STRG.33537.1 loc:chrX|125404946-125405412 - Another Rps12 pseudogene.     
STRG.20579.1 - Some Pfam Atpase coverage (in bacteria), predicted TM domain... 11x coverage - modest conservation across gene structure.     






