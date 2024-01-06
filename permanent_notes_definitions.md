# Document Outline
Herein will be permanent notes, musings and writings regarding the analysis of
the datasets contained within this folder. It will include definitions,
thoughts, and any other ideas which are necessary for the reproducibility, and
understanding, of the current data.    

- For further completeness, unknown transcripts were cross-referenced against
the RefSeq annotation also. Any overlap at all (even 1 base) indicates
transcriptional coverage that has already been detected.    

- DeSanta (2010) enhancer paper uses 20kb upstream as a cut off for a distal
enhancer. I am taking the same approach for the transcripts which don't have a
CRE signature, even though some enhancer may be located further away. I would
rather be conservative and careful.    


**UTR's and Polymerase run-on**
- Definition of a 3'UTR extension and/or polymerase run-on is to be within 5kb downstream of another gene, being transcribed in the same direction
- Definition of a 5' UTR extension is to be within 5kb upstream of another gene, being transcribed in the same direction.
- Transcripts which appear to be UTR extensions, and polymerase run on, but contain poly-A tails... what to make of these? Internal priming sites?

**Some Definitions below**         
* 3' UTR extension = 5KB downstream on SAME strand
* 5' UTR extension = 5KB upstream of a gene on the SAME strand - this will need
to be carefully considered as it may also be a Promoter associated
transcript/PROMPT and so CRE-promoter layering is beneficial, CPG islands.
* Upstream of an adjacent gene by 20kb on the same strand is classified as an
enhancer , put transcripts further than 20kb which contain a CRE-enhancer mark
are also included in this set.
* Pseudogenes are those which have high scoring BLASTN hits to refseq-known -
* Intergenic nc transcripts which are probably described as lincRNAs - TUFs/Transcriptional Units.


Minimap version used for all analysis - 2.17 

