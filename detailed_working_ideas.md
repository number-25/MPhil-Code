# Document Outline
This document was born from the need to catalogue my ideas and strategies for
analyzing the remainder of the directRNA dataset. I was posting everything to
the slack message board, until it because overwhelming and confusing. 

* In silico fasta read fragmentation - Take pseudogene co-ordinates, as well
  as their extant genes, digest them with a short read simulator, then map with
  short mappers, at various depths, and see if they can be detected - do the
  same with a long-read simulator.      https://github.com/liyu95/DeepSimulator
  https://onestopdataanalysis.com/read-simulator/
  https://github.com/schmeing/ReSeq https://github.com/HadrienG/InSilicoSeq
  https://github.com/wanyuac/readSimulator https://github.com/qasimyu/simuscop
  https://www.biostars.org/p/1852/ 
  GemSIM: general, error-model based simulator of next-generation sequencing data
  https://academic.oup.com/bioinformatics/article/31/17/2778/183245     
  http://alumni.cs.ucr.edu/~liw/rnaseqreadsimulator.html
  https://www.biostars.org/p/288294/     
  https://confluence.sammeth.net/display/SIM/4.1+-+Gene+Expression    
https://academic.oup.com/nar/article/40/20/10073/2414449
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3450-9

* Internal priming detection algorithm - integrate it into the basecaller, read until?    
* DeepTools has some pretty interesting programs for plotting chip-peaks, quality control computes.     
* If short-reads can in fact detect these - then perhaps there are more optimal
fragmentation strategies which favour certain hexamer biases for mapping? Or is
it a function of read depth - they can be detected, but must be sequenced
extra-deep to achieve the degree of resolution needed?    


## Article links and readings    
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2766791/ (Hidden Markov models intro)     
- NCBI eutils cookbook - need to get more comfortable with this https://github.com/NCBI-Hackathons/EDirectCookbook    
- Tutorial: Python Regex (Regular Expressions) for Data Scientists - https://www.dataquest.io/blog/regular-expressions-data-scientists/   
- Evolutionarily conserved elements in vertebrate, insect, worm, and yeast genomes   
- Notion - Notes/NotesBetween Computers - for RNA-seq and visulisation tools 
- Read proudfoot poly-a tail review "Ending the message: poly(A) signals then and now" when dealing with poly-a data.
- Poly-A signal plotting examples. "Precise gene models using long-read sequencing reveal a unique poly(A) signal in Giardia lamblia"

## Dave Tang    
Defining genomic regions (how to create intergenic, intronic bed files etc.)    
https://davetang.org/muse/2013/01/18/defining-genomic-regions/     
https://davetang.org/muse/2012/08/07/sequence-conservation-in-vertebrates/
Sequence logos https://davetang.org/muse/2013/01/30/sequence-logos-with-r/
PWM https://davetang.org/muse/2013/10/01/position-weight-matrix/

## Evolutionary Analysis 
- Retrieving Phylop Scores https://www.biostars.org/p/86847/
- Conservation score for a piece of human genome (PhastCons and PhyloP) https://www.biostars.org/p/152454/#152468    
- getting a conservation score from UCSC if I have exon coordinates? - https://www.biostars.org/p/242484/#242486   
- Article - "Evolutionarily conserved elements in vertebrate, insect, worm, and yeast genomes"    
- Detection of nonneutral substitution rates on mammalian phylogenies - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2798823/   
- Of those that may be conserved and in close proximity to protein coding
genes- how many are hitchhiking, in for the  ride, and haven’t been recombined,
are in fact still linked.     
- align orthologues, then compare to neighbouring neutral sequence (ancestral
repeats) of similar GC composition - to see if more or less pronounced
evolutionary rate.     
- Are the regions of conservation overlapping with the TF binding site.      
- RNAcode on high confidence protein coding candidates - create multiple sequence alignment, then run RNAcode on it, in additional to PhyloCSF    
- Conservation around splice sites and promoter region    
- Perhaps average phyloP score across gene body? Compared to protein coding genes - but how to highlights conserved TFBS/regulatory binding sites?    
- Sequence conservation in vertebrates https://davetang.org/muse/2012/08/07/sequence-conservation-in-vertebrates/
- Perhaps average phyloP score across gene body? Avg. across body... compare to protein coding genes avg. - Also plot phyloP score across gene body to posit the notion that the 3'UTR typically has elevated conservation signal.,
- Use putative PC, NC genes, and AR as outgroups. 


## Protein Coding set
- Map to proteo-genomics datasets - look for assistance in APRIS, tress work - SQANTI information (pride database) for unique discriminatory peptides.
- I fear that this will be my greatest difficulty - attempting to parse these
and understand what to keep and what to discard. Expression measures will
have to inevitably inform my decisions. 


## Coverage plots and other R plots
- Make an expression graph of the non-coding set - have the ability to color
code the enhancers to show that most of them have such minimal expression
levels     


## Poly-A tail and Poly-A site
- Poly-A stretch in their sequence, rather than having a poly-A tail? Also, is the polyA site within the last 30nt of the RNA?
- Internal priming of A rich sequences, thus picking up aretefactual reads.    
- Enrichment of PAS in enhancer associated RNAs, and even within the junk transcripts.     


## Enhancer RNAs  
- DeSanta 2010 used 10kb from 3' end for transcriptional “run off” filtering.
Also “First, we assigned predicted transcribed enhancers to adjacent coding
genes if distant from them less than 20 kb.”    
- Using enhancer marks as positive evidence for enhancer associated RNAs - H3K4me1, H3K27ac, p300, CBP binding 
- if something is monoexonic, low expression, and in very close proximity to
another genetic element, chances are, it is an unproductive transcript or
simply, a 'transcriptional unit
- Ser5P  + Ser2P profiles in mice - if any datasets exist, they would be super helpful. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73444 
https://link.springer.com/article/10.1186/1471-2164-12-516 - net-seq mouse cells RNA Polymerase II Phosphorylated on CTD Serine 5 Interacts with the Spliceosome during Co-transcriptional Splicing
- GSE43390 and GSE27037 GEO accesion numbers for GRO-seq - looking for pol-II pausing.     
- PWM for TFs on the transcripts believed to be enhancers, promoters etc.
- If a transcript has a PLS, and also a pELS, it is most likely a PROMPT,
especially given it runs opposite to the annotated gene, as contains a CPG
island.      
- JASPAR database + PWM MEME suites.     

## Promoter associated transcripts 
- Containing CpG island, refTSS signatures, high H3K4me3 
- Independent, intergenic promoter platforms 
- Similar to annotated cleavage sites at TESs of genes, these upstream
antisense cleavage sites are associated with the PAS located at the expected
position, about 22 nucleotides upstream of the cleavage site” - PROMPTS in my
dataset, the ones with poly-A sites, and poly-A tails. From “Promoter
directionarlity is controled by U1…”
- It appears that 3kb upstream and divergent (antisense) appears to be the cap on the production of these transcripts. 
- Is there a sweet spot in the distance from the closet gene, and the potential for it being designated a prompts? e.g. after 5kb do they then begin losing the PAS enrichment etc.? Are there particular discirimators to aid in the annotation process.     


## Splicing and Exon content
- How many are spliced and multi-exonic? Do they product alternative transcripts? 

## 3' UTR extentions and polymerase run on/rippling
- Transcription termination signal/motifs to identify separate transcripts from UTR extensions    


## Position Weight Matrix, Motif analysis and MSA presentation
- https://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf 
- https://cran.r-project.org/web/packages/ggmsa/vignettes/ggmsa.html
- https://hartleys.github.io/QoRTs/archive/v1.0.7/Rhtml/makePlot.cigarOp.byCycle.html
- WebLogo http://weblogo.berkeley.edu/logo.cgi  
- https://davetang.org/muse/2013/01/30/sequence-logos-with-r/ 
- https://davetang.org/muse/2013/10/01/position-weight-matrix/
- https://www.ebi.ac.uk/Tools/msa/mview/    

## Pseudogenes   
- For the hypothesis that TE's are using the 3' UTR of existing genes to
facilitate their transposition (by perhaps using PAS signals to avoid
degradation) - then I must find a lineage in which the orthologue is flanked by
the same TEs (?), but hasn't been mobilized yet, but then, upon jumping, mobilizes the UTR.     
- Another big question is, given only a small fraction of the 3', last exon is
being detected as pseudogenic, how does the gene get transcribed ? Need to look
for pol-2 sites, promoters etc. 5' -> 3' degradation ? the idea that we're catching transcripts JUST before they are degraded?     
- Another hypothesis is that the 3'UTR of these bonafide gene is in fact viral in origin?  
- I have noticed that some has a repeat masked string of poly-As - could this
in fact be donated by the TE upon insertion?    
- Are the 3'UTR protecting the TEs from degradation - hijacking this. 3'UTR
have anti-degradatio signals within them allowing the pseudogenes to migrate?       
- GO search of the genes which have had their 3' ends copied. Anything
particular about the sequence content which may make it more likely for these
"types" to be transposed?      
- To show that they have transposed - you'd have to find an orthologue which borders the same TE in another species. 
- Many of these intersect with non-mouse mRNAs from ncbi, indicating they are active in other animals?     
- Why are the 3'UTRs detectable? : 3'UTR is region with some of the highest
conservation signature, meaning its attrition will be longer than the rest of
the sequence. One could intersect these wit the LOD conservation blocks as most appear to have them.       


## Structural RNA and small RNA    
- Infernal scans, RFAM
- http://eddylab.org/software/snoscan/snoscan.README http://bio.lundberg.gu.se/srpscan/rnascan_doc.html tRNA scan lowe lab 
- For predicted tRNAs, those with low scores (under 30), if they overlap or are
in very close proximity to B2 sine elements, it may be a positive hit due to
the native tRNA sequence within such elements - so this is needed for
confirmation.      

## Repeat associated and Paralogues / High Copy Number
- Are genes with high copy number and elevated repeat content (masked lower
case) more likely to be unn-annotated - is it a function of this?     

## LincRNAs
- Canonical splice junctions?   
- Use H3K36me3 from Bing Ren (UCSC) to overlap with potential lncRNAs https://www.biostars.org/p/150036/     
**to merge all H3K36me3 marks cat tmp.reftss tmp.nc tmp.pc | bedtools merge -i - > EPD-REFTSS-merge.bed**     


## Bambu
- https://stackoverflow.com/questions/59792855/loading-fasta-file-in-r-faster-than-when-using-read-fasta-from-seqinr
- https://github.com/GoekeLab/bambu
- https://rdrr.io/bioc/Rsamtools/man/scanBam.html
- https://kasperdanielhansen.github.io/genbioconductor/html/Rsamtools.html
- https://www.biostars.org/p/394233/ 


