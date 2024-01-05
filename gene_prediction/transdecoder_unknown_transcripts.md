# Document Outline 
The possibility in which 'unknown' transcripts within our Direct RNA dataset
may contain ORFs, ought to be investigated. Here, this attempt is embarked on.

# Aims
* Use TransDecoder for the first time.  
* Dedicate some time into thinking about what other strategies may be employed
to investigate this notion of ORFs within unknown transcripts.    
* Use InterProScan for the first time.   
* Use Trinotate for the first time.  


# Analysis    

## TransDecoder
($FA is the .fasta sequence of the unknown transcripts).  
`TransDecoder.LongOrfs -t $FA`  

Download UniRef90 from UniProt, and the latest Pfam HMM, in order to undertaken homology search, the
results of which will be integrated into the final predictions.  
```
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz  
ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz  
```   
**Redid TransDecoder analysis after I removed pseudogenes, tRNAs, and transcripts arising from ambigious chromosomal scaffolds.**    

Run the Blastp search against UniRef90.  
**Need to create a blast-db - most likely in the cluster.**   

Transfer over uniref90.fasta files to cluster.  
`rsync -av -P
/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Proteomes/UniProt/uniref90.fasta
uqdbasic@tinaroo.rcc.uq.edu.au:/30days/uqdbasic/input`  

Create the blast-database for uniref90.fasta.  
```  
#!/bin/bash   
#PBS -A UQ-QBI   
#PBS -l select=4:ncpus=3:mem=70GB   
#PBS -l walltime=36:00:00   
#PBS -N makeblastdb  

module load blast2/2.10.1   
INP=/30days/uqdbasic/input/uniref90.fasta   

makeblastdb -in $INP -input_type fasta -dbtype prot -title uniref90
-max_file_sz 3GB -parse_seqids -out /30days/uqdbasic/output/uniref90-blastdb   
```   
PEP=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Anno
tation/Trinotate/TransDecoder/RC-and-EXT/true-unknown.fa.transdecoder.pep
Downloaded locally from cluster.  
Run local blastp search using the newly built uniref90 database.    
```
PEP=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/TransDecoder/RC-and-EXT/unknown-txs.fa.transdecoder_dir/longest_orfs.pep   

# Be very, very particular and observant of the database name - the .fasta, and
the resulting .x database files must all have the same prefix, and blast should
be pointed specifically to this prefix, which should be identical between files
e.g. 'uniref90-blastdb.*'

UNIREF=/media/labpc/Disk-2/Bioinformatics-Computational/blastdb/Uniref90-db/testing/uniref90-blastdb

blastp -db $UNIREF -query $PEP -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > /media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/TransDecoder/RC-and-EXT/blastp.outfmt6
```

Run hmmscan (hmmer).   
```
HMM=/media/labpc/Disk-2/Bioinformatics-Computational/Reference-Proteomes/Pfam/Pfam-A.hmm  
hmmsearch --cpu 6 --domtblout pfma.domtblout $HMM $PEP   
```  

* In order to understand the output of hmmscan with Pfam, one must read the user manual written by Sead Eddy.  

`~/Desktop/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t $FA
--retain_pfam_hits pfma.domtblout --retain_blastp_hits blastp.outfmt6`  

* **Of the newly predicted ORFs/proteins, which ones overlap with existing
enhancers, and transposable elements, and may thus already be accounted for.
In particular the transposase and reverse transcriptase activity elements of TE's will surely be detected in
these searches.** 

## InterProScan 
Interproscan is one of the automated gene annotation pipelines utilised by
EMBL-EBI. I have not used it before, however, it shows promise, especially for
protein coding genes. The goal here is it process the unknown.fasta files
through it, then compare the output with that of
trinotate/transdecoder/pfam-blastp-hmmer, and parse the output into confidence,
likelihood categories.  

interproscan.sh -d $OUT -gotrems -i $IN-path -t n 

```  
IN=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/true-unknown.fa 

OUT=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/InterProScan/

./interproscan.sh -d $OUT -goterms -i $IN -t n 
```   

The pipeline ran correctly, although I had to disable tmhmm search.  

Sort by highest e-value - create a separate file with basic information in it.   
`cat *.tsv | cut -f1,3,4,5,9 | awk '$5!="-" {print $0}' | sort -k5,5 -g > true-unknown-sorted-top-hits.lst`   


## Trinotate 
In an attempt to utilize trinotate to annotate the unknown transcripts that we
have, I have to create a custom gene-to-transcript file, as this is a necessary
file for input. Typically, one derives this file by running the entire trinity
pipeline. As I am take a different route, I must supply this file my self.
Stress less, as it appears it can be created from the output files from
TransDecoder, as described
(here)[https://github.com/Trinotate/Trinotate.github.io/issues/7].  

Manipulate the .gff3 file to extract the necessary columns described above. The
gene name, and the transcript name given by TransDecoder. Gene-to-transcript.map.    
```
cat true-unknown.fa.transdecoder.gff3 | awk '$3=="gene" {print $0}' | cut -f 9 | awk -F"~~" '{print $1}' | sed 's/ID=GENE.//g' > trino.txt

cat true-unknown.fa.transdecoder.gff3 | awk '$3=="gene" {print $0}' | cut -f 9 | awk q:-F"~~" '{print $2}' | sed 's/;.*//' > trino2.txt

paste -d"\t" trino.txt trino2.txt > gene-to-transcript.map
```   

Now that all of the recommended [modules](https://github.com/Trinotate/Trinotate.github.io/wiki/Software-installation-and-data-required) are installed, I will proceed to run the data through the programs. 

Blastx.   
```
makeblastdb -in uniprot_sprot.pep -dbtype prot -out /media/labpc/Disk-2/Bioinformatics-Computational/blastdb/uniprot-sprot-db/   

blastx -query $FA -db $UNIPRO -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6    
```

Hmmscan with modified output.   
```  
PEP=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Anno
tation/Trinotate/TransDecoder/RC-and-EXT/true-unknown.fa.transdecoder.pep   

hmmpress Pfam-A.hmm     

hmmscan --cpu 8 --domtblout TrinotatePFAM.out $HMM $PEP > pfam.log   
```   

tmhmm scan. (Have to call tmhmm in it's absolute path).      
`/home/labpc/Bioinformatics-Programs/tmhmm-2.0c/bin/tmhmm -short < $PEP >
/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/tmhmm.out`    

Run rnammer to search for rRNA genes.   
`/home/labpc/Bioinformatics-Programs/Trinotate-Trinotate-v3.2.2/util/rnammer_support/RnammerTranscriptome.pl
--transcriptome
/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/true-unknown.fa
--path_to_rnammer /home/labpc/Bioinformatics-Programs/rnammer/rnammer`    

* rRNA was not detected, .gff is empty.   

Load data files into Trinotate SQL database.   
```   
admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate  

MAP=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/TransDecoder/RC-and-EXT/gene_trans_map    

FA=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/StringTie/long-
read/RC-and-EXT/gffcompare/unknown/no-annotation-truly/true-unknown.fa   

PEP=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Anno
tation/Trinotate/TransDecoder/RC-and-EXT/true-unknown.fa.transdecoder.pep    

BLASTP=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-A
nnotation/Trinotate/TransDecoder/RC-and-EXT/blastp.outfmt6    

# Remove gene ID's from the blastx.out which were not retained in the transdecoder output, as this is causing confusion in the trinotate sql loading.  

cat gene_trans_map | cut -f1 > tmp.lst   
fgrep -w -f tmp.lst ../../blastx.outfmt6 > ../../mod-blastx.outfmt6    

BLASTX=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/mod-blastx.outfmt6   

PFAM=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Ann
otation/Trinotate/TrinotatePFAM.out   

TM=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annot
ation/Trinotate/tmhmm.out    

SP=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/signalp.out   

./Trinotate ./admin/Trinotate.sqlite LOAD_swissprot_blastp $BLASTP LOAD_swissprot_blastx $BLASTX LOAD_pfam $PFAM LOAD_tmhmm $TM LOAD_signalp $SP

./Trinotate ./admin/Trinotate.sqlite report -E 1e-3 > /media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/
Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Annotation/Trinotate/trinotate_annotation_report.xls
```   

# Some top-scoring examples extracted from both pfam and blastp outputs.   
STRG.675.1    


## High scoring protein coding candidates 
STRG.19052.1 - Cov: ~68, InterPro: PANTHER PTHR11931-PHOSPHOGLYCERATE
MUTASE-10:262-3.3E-174, Pfam: PF00300.24 with 4.1e-36, UniProt: A0A4X2LNU3 -
98% identity - eScore of 0.0, Positionally conserved in many species,
vertebretes, etc. in UCSC. Has a very, very similar pseudogene (pgam) on chr7,
so be weary.   

STRG.3804.1 - Cov: ~52, InterPro: Gene3D - Classic Zinc Finger - 54:93 -
4.0E-14, Pfam: PF13465.8 - 4.2e-08, Uniprot: Q8BYD9 - 56.522 - 5.17e-25 - May
occur in multiple domain repeats, Not super conserved, has some TE's within it
- overall looks ambiguous in BLAT, First reading frame 5'-3' has substantial
domain architecture as shown in SMART - homology is TEC gene ENSG00000135605,
second and third 5' - 3' reading frame has single classic Zinc finger
Domains. Has a poly-A cluster within it.  

STRG.16252.1 - SpoIIP sporulation protein.    


STRG.25430.1 - UPI0007ECD13E - Many repeats within it (MANY!!), has a predicted
poly-A cluster (weak), and strong ATAC peaks.    

STRG.29286.6 and STRG.29286.4 - multi-domain protein? a remnant? Pseudogene? fascinating.  

STRG.22190.1 - Decent PhyloP - decent expression level (~40), has a poly-A
cluster within it, uniref hit, has a predicted signalP hit.  

STRG.21070.1 - Complex locus - many repeats - many TE influences - decent
expression ~36 - high multiple transdecoder hits, multiple uniref hits.    

iSTRG.29298.1 - FAM pfam domain (spermatogenesis) - internal repeat domain has a
BLAST hit to a cell surface glycoprotein - has a protein architecture in SMART
(two repeats and TM domain), has predicted transmembrane regions in tmhmm.out.   

STRG.29314.1 - Repeat associated? Pseudogenised? Complex region resolution -
Has many many hits, cell surface glycoprotein? Use this as an example of
complex region. Has several, very similar repetitive elements surrounding it.
Interpro-scan!!  

Something such as STRG.7033.1 is bizarre... upstream of a gene.. has some
methylation activity, PEPcase pfam homology, but PEPcase is not found in humans
or fungi?? Uniref hits, predicted signal peptide, predicted transmembrane
domains... predicted ORF... 12 coverage score...    

STRG.10850.1 - Fibroblast growth factor Orthologue? Predicted protein in
chlamydia, weak ORF? Clearly maps to cente of a LINE element? What the fuck?
Align the fibroblast ortho to the line element? 


#### Potential paralogues 

STRG.8996.1 - Gm21103 gene copy or pseudogene - has many hits across the chromosome 14, potentially indicating multiple repeats and duplicates.   

STRG.3249.2 - Potentially a Gm6729 gene copy or pseudogene - has many hits
across the chromosome like aboe. Has a hit to the intron of the parent gene
which is being strongly selected against and purified.   

STRG.18562.1 - A gene copy/pseudogene of Cid3? Covers almost the entire
transcript in the parent gene - not sure about it's gene structure at the novel
locus.    




#### Protein coding pseudogene candidates.

STRG.30314.1 - Arl5a pseudogene, 3' behaviour - Exact same TE elements within sequence.  

STRG.6583.1 - 2900097C17Rik pseudogene as above, flanked by TEs in novel locus, multiple hits. 

STRG.12738.1 - Another Ppp1r2 pseudogene which hasn't been annotated yet.    

STRG.33537.1 - Another annotation for Rps12 pseudogene which hasn't been annotated yet.     

STRG.17533.1 - Nemf pseudogene? 3' behaviour.    

STRG.26035.1 - Akirin1 3' end homology (1kb homology, primary transcript is 5kb)    

STRG.16503.1 - Dusp7 pseudogene as above, novel locus flanked by TEs    

STRG.21003.1 - Eif4b pseudogene as above - multiple hits across genome (multiple duplication/transposition events.) Flanked by TE     

STRG.21632.1 - Zfand5 pseudogene with 3' end homology - (1kb homology, primary transcript is 5kb) - large part of it is repeat masked. 



#### Non-coding 

STRG.3575.3 - good example of a very lowly expressed, multi exonic transcript,
which has phyloP and phastcons signal, cis-regulatory marks associated with it
- but such, such low expression. What to make of these?    

STRG.9115.1 - another pseudogenic example - very high homology to Csde1. Gooodo 

STRG.6966.1 - Example of a very lowly expressed transcript, which overlaps
ENCODE cis-regulatory elements as a distal enhancer. Has substational
conservation 

STRG.25631.1 - Another (high expression albeit) pseudogene  **reavaluate this**    

STRG.19405.1 - An example of a likely putative enhancer RNA 

STRG.17817.1 - Another example of a very similar pseudogene - to Gosr2 - very high homology to the 3' end of it.  

STRG.20865.1 - Has the highest coverage in our dataset (by a looong margin).
Extensive and deep conservation (Phylop), has been detected before in mouse
cDNA screenings (riken etc.), has RNAcentral sequence similarity to a rRNA, has
BLASTn hits to rat matrix metallopeptidase 16 (Mmp16) and ~50% alignment
coverage. A deactivate gene? Protein coding turned non-coding? Frame
translation gives many transmembrane domains inside it. Compare to Mmp16, see
where the start and stop codons are, and perhaps observe whether frameshifts
and PTC were introduced, start codons destroyed. Other non-mouse mRNAs detected
- mostly from species where the brain was the tissue being sequenced...       

STRG.3789.1 - Contains a miRNA within it - mir-684, also matches the intronic
sequence of another gene, giving cause to believe that this may be another
pseudogene - Dusp19 or Smt3h2-ps4.     


STRG.12739.1 - Tulp4 Pseudogene, has similarity to the 3' end of the last exon.
- As above. Conservation signal which is maintained in other loci - upstream of
pseudogene in both novel locations.      

STRG.12732.1 - Tulp4 Pseudogene, or lncRNA? (Two pseudogenes found) - As above.
It's other two locations both also have very similar high conservation signal.     

STRG.3955.1/2 - Spin1 Pesudogene, once again, to the 3' end of the last exon - As above.

STRG.21206.1 - Strn Pseudogene - small section of last exon.        

STRG.9525.1 - Fabp5- pseudogene - Bizarre, also has high identity to the
introns of other genes - all around 180-200bp in length, and maps quite highly
to other genomic regions which no annotations. **Many of these are flanked on
one side by simply-repeats??** - The original match has a low-mid conservation
signal which all surrounding regions do NOT have. Is this a structural element? Infernal!  

STRG.9001.1 - Very strange - has many, many hits across chr14, including hits
to the exons of protein coding genes - an undocumented transposable element?
What could it be?     

STRG.31179.1 - Lrrc58 Pseudogene - Once again, to the 3' end of the last exon - As above. May also have another pseudogene on chr2. 

STRG.17463.1 - Has conservation signal and similar mRNAs in other species, has
an exact duplicate which maps downstream to the intron of a pseudogene, and
also maps to the 3' end of another gene Ptprj, which is located downstream on
the same chromosome. Spooky.      

STRG.17042.1 - Investigate further.   

STRG.32708.1 - Very bizarre, the region across it's intron has no information -
no GC content, no chip, conservation, no atac, no SNP insight, NOTHING across
the region in UCSC. On it's 3' end, it has two CPG island, followed by strong
promoter signals - indicating, it could be serving as an enhancer, leading into
a promoter region? Has high coverage.     

STRG.7109.1 - Clearly has enhancer like qualities - same strand from gene, CRE overlap, mild conservation.    

STRG.7196.1 - A notably bizarre example - Has over 100x coverage, located in
very open ATAC regions, has homology to an un-annotated element in the
mitochrondrial region, has conservation, and where it gets wierd, is that the
other highly similar regions of the genome (via BLAT) also all have open
chromatin loops and conservation signatures similar to pseudogenes, and are
neighboured by LINE TEs, but the region of ATAC has been purified of TE sequence? 
Structural element?     

STRG.25631.1 - Super high scoring, (long) hit pseudogene for Rprd1a - Once
again, homology to the last exon!!! Virtually no coverage.   

STRG.28365.1 - As above, simply for a different gene - Ero1

STRG.19121.1 - Nudt21 Pseudogene, just as above!   

STRG.18827.1 - A nice example of a multiple duplication event in which a short
~500nt element is located multiple times across a speficic region of a
chromosome - this eludes the status of functionality and makes obvious the
pitfall of designating functionality based only on transcription.    

STRG.32048.1 - Mthfs Pseudogene - Once again, on its last exon,linked to the poly-A tail.    

STRG.8738.1 - Potentially an intergenic transcript - has CpG site, and is 8kb away from neibhouring gene.. hmm 

STRG.28674.1 - Zbt18 Pseudogene - Same relationship as the others - last exon,
linked to poly-A tail. Has potentially other pseudogenes in genome, one within
an intron of an existing gene - snoRNA? smallRNA?     

STRG.10193.1 - Nampt pseudogene - As above - 3' end,     

STRG.10697.1 - Nadk2 pseudogene - as above....     

STRG.1244.1 - Fam119x pseudogene - as above... again     

STRG.6437.1 - Nice example of a duplication event (it seems) which results in a
palindrome or mirrored like sequence which contains virtually identical phylo
marks, repeat elements, etc. - perhaps this is the incompleteness of the genome
assembly which has led to this? It's mirror is farther away (screenshots are on
the desktop).    

STRG.1758.1 - Dok7 pseudogene/homology? This time on the 3' end, but not completely towards the end..?     

STRG.4019.1 - Perhaps a pseudogene? Maps upsteam of Actr2, yet has pretty high homology to the 3' end of Kpnb1 (how to resolve this?!)     

STRG.21297.1 - Vmp1 pseudogene - follows same pattern as above 3' end.     

STRG.22552.1 - Hnrph2 pseudogene - as above.    

STRG.33513 - Foxj3 - as above! Whatttt?????   

STRG.19335.1 - U2af2 pseudogene - as above.     

STRG.1189.1 - Ap3s1 pseudogene - on 3' last exon, in the middle of it. Also
found on ChrX, Chr8 once intergencally, and again withn the intron of Mcph1 -
are these perhaps somesort of small RNAs? structural? I am indeed lost for
words....    

STRG.14624.1 - Shoc2 pseudogene, as above.     

STRG.6262.1 - Found within the intron of Invs - flanked by a LINE (transposed
with it?) and pseudogene on other end. Same pattern in Gm47882. Also has high
homology to the 3' end of the last exon of Hnrpf. Bordering a pseudogene in 3 
other genomic loci also - isn't it strange that it's neighbouring them yet they
are not related?    

STRG.32988.1 - Fem1b pseudogene matching the same 3' last exon patter. Flanked by LINES at its novel locus? A coincides? 

STRG.3479.1 - Rfx7 pseudogene as above.     

STRG.14676.1 - Cctn pseudogene as above, also has homology (95%) to another intergenic locus, 

STRG.17674.1 - Kif5b psuedogene, as above. Aaaaand wtf! 95% homology to an intronic region of Cplane1 and flanked by TEs!     

STRG.9208.1 - Kif3c pseudogene, as above. 

STRG.12346.1 - Azi2 pseudogene, as above.    

STRG.453.1 - Eif3a pseudogene. Also has 95% homology to the intron of Pcnp.    

STRG.1127.1 - Caprin1 pseudogene as above. Flanked by TEs in novel locus, and at a lower homology locus it is right next to a tRNA gene.    

STRG.8189.1 - Prrc2c pseudogene as above, multiple hits across the genome - structure again? transposition?    

STRG.33541.1 - Hnrnph2 pseudogene as above.flanked by TE at novel locus. Has has high homology hits elsewhere which are flanked by TEs.   

STRG.26398.1 - March6 pseudogene as above - flanked by TE at novel locus.

STRG.25769.1 -  Cd2ap pseudogene - flanked by TEs at nove locus.

STRG.31221.1 - Lrcc58 pseudogene - as above 3' end and last exon. Multiple high homology hits across genome.   

STRG.12737.1 - Another potential pseudogene/duplication - Has homology to the
first intron and exon of Dynlta. In this region, it overlaps a CpG island, in
it's novel locus 1 it has two CpG island very close with positive promoter
signaturs. In its second two loci, it borders a copy of the same pseudogene, as
well as it's close proximity to CpG islands - this duplication like signature.    

STRG.11364.1 - Novel locus is flanked by LTRs / TEs, second novel locus is flanked by LINE,
3rd by LTR, 4th by LINE, homology to intron of Cpne4 - flanked by decaying LINE. Multiple hits across genome - all transposable?        

STRG.14810.1 - C33007P06Rik pseudogene - follows same patterns are others. Flanked by small SINE on one end.     

STRG.3981.1 - Hmgb1 pseudogene (flanked by small TE) - as above. Multiple hits across genome. Flanked by Hmgb1-pseudogene in other loci.    

STRG.14091.1 - As above. 

STRG.19880.1 - Rapgef2 pseudogene as above - Several homologous locations - all
individual insertions at each site are flanked by identical upstream TE.. 

STRG.26577.6 and .1, .5, .2, .4 - Actg1 pseudogene - As above. Has a hit within intron
of other genes. Multiple hits spread across the genome - flanked by HUGE TE in
an observed locus.        

STRG.33795.1 - Pgrmc1 pseudogene - As above.     

STRG.14489.1 - Set pseudogene as above. Multiple hits across genome.      

STRG.30854.1 - As above. Novel loci are flanked by TEs.     

STRG.29285.1 - Larp1 pseudogene - as above. Muktiple hits across genome.    

STRG.20256.1 - Morf412 pseudogene - as above. Flanked by TEs in novel loci -
flanked by similar pseudogenes at most locations... Are these all "different"
pseudogenes, or the same ones, in various states of decay.     

STRG.20901.1 - Foxk2 pseudogene - novel locus flanked by repeats.    

STRG.33756.1 - Tmtc3 pseudogene as above -  Flanked by huge TE. 

STRG.31971.1 - Likely a strongly decaying pseudogenic insert - has homology to some introns in bonafide genes.  

STRG.15111.1 - Stgal5 pseudogene with 3', last exon signature. Flanked by TEs.   

STRG.3955.2 - Spin1 pseudogene as above - one long 1.5kb homologue, with multiple smaller loci - most seem flanked with TEs.     

STRG.19464.1 - Nice example of filling out an existing pseudogene model - has homology to last exon of Eef1a1. Multiple hits across genome.     

STRG.32950.1 - 1810037l17Rik pseudogene - begining of last exon - has 4 novel
hits across genome, one right next to annotated pseudogene (the rest of it?).   

STRG.25001.1 - Homology to the 5' and 3' heads of Khsrp and its ajdacent gene - at novel locus it borders a pseudogene.  

#### Those to omit - 
STRG.27809.1 - already a known pseudogene. 
STRG.21125.17 - TE element
STRG.16551.1 - known pseudogene
STRG.9426.1 - TE element 
STRG.24664.1 - PROMPT     
STRG.8588.1 - TE
STRG.17042.1 - TE or intron remnant... not much to say
STRG.21065.1 - Already annotated pseudogene - small homology
STRG.27212.1 - TE   
STRG.19177.1 - TE 
STRG.23624.1 - TE
STRG.17334.1 - PROMPT
STRG.16015.1 - Enhancer/Regulatory   
STRG.16535.1 - TE   
STRG.11807.1 - PROMPT + TE    
STRG.20185.2 - Random TUF   
STRG.26848.1 - TE    
STRG.31987.1 - TE 
STRG.28543.1 - TE
STRG.21275.1 - TE
STRG.16936.1 - TE
STRG.204.1 - TE
STRG.10908.2 - TE    
STRG.19904.1 - TE  
STRG.26684.1 - SINE   
STRG.13827.1 - TE LTR
STRG.19574.1 - TE    
STRG.29962.1 - TE    
STRG.8283.1 - TE  
STRG.30837.1 - TE  
STRG.25352.1 - TE 
STRG.23624.1 - TE 
STRG.26149.1 - TE   
STRG..7176.1 - TE  
STRG.13778.1 - TE   
STRG.28954.1 - TE   
STRG.6244.1 - TE    
STRG.1753.1 - TE  
STRG.9182.1- TE   
STRG.17182.1 - TE 
STRG.17337.1 - TE 
STRG.20846.1 - TE   
STRG.23083.1 - TE 
STRG.205.1 - TE 

