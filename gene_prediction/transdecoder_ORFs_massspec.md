# Document Outline     
Herein, I will entertain whether the predicted open reading frames from our
dRNA data, have high confidence mappings to the mass-spec data that Wei-Siang
has generated, attempting to establish some preliminary links between these
predicted peptides and experimentally derived ones.    

## Aims
* Map the predicted ORF .peptides to the Mass Spec derived peptides for the first time with BLASTP.     
* Parse an excel document and format into a fasta file.     


## Analysis 

### Excel file parsing.    

Convert the excel format files (xlsc?) to .fasta format with the sequence name in the header and the peptide sequence in the body.      
One can use awk for this, as pointed out (here)[https://www.biostars.org/p/423573/].    

The sequence name is in the 7th column, and the peptide sequence is in the 13th.       
Save the excel file as a .csv. 

Convert from a comma separated to a tab separated file. Cut out only the
relevant columns. Remove the header. Print the columns using a simply awk
script.    
Some of the columns are not formatted correctly (classic), so I will also need to extract the peptide sequences from the 14th column.   

```    
less -S Nuclear-EXT-replicate-1-peptide-summary.csv | awk -F, '$13 ~ "[A-Z]" {print ">"$7"\n"$13}'    

touch tmp.tmp 

less -S *.csv | sed '1d' | awk -F, '$14 ~ "[A-Z]{3,}" {print ">"$7"\n"$14}' >> tmp.tmp     

less -S *.csv | sed '1d' | awk -F, '$13 ~ "[A-Z]" {print ">"$7"\n"$13}' >> tmp.tmp     
```    


Clean poorly formatted .fasta headers, remove duplicate entries, and add unique
identifiers to duplicate fasta accesions (with different sequences).       
`less tmp.tmp | sed "s/\;.*$//" | seqkit rmdup -s - | awk '(/^>/ && s[$0]++){$0=$0"_"s[$0]}1;' - > Mass-Spec_WS-Nuclear-RC_EXT-pep.fa `     
* 14193 peptide sequences present.        

Old attempt which didn't work correctly - errors showed up when building the blast database - duplicates and wrong values.    
`less -S Nuclear-RC-replicate-1-peptide-summary.csv | sed 's/,/\t/g' | cut
-f7,13 | sed '1d' | awk '{print ">"$1"\n"$2}' >>
Mass-Spec_WS-Nuclear-RC_EXT-pep.fasta`  


### Create a blast database and search with blastp.  

Add unique identifiers to each peptide .fasta record. From
(here)[https://bioinformatics.stackexchange.com/questions/979/how-to-append-numbers-only-on-duplicates-sequence-names]    
`awk -iinplace '(/^>/ && s[$0]++){$0=$0"_"s[$0]}1;' file.fa`    


```
makeblastdb -in ../mass-spec-outfiles/Mass-Spec_WS-Nuclear-RC_EXT-pep.fa -input_type fasta -dbtype prot -title WS-MassSpec-pep -out ./    

PEP=/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/Gene-Pred-Anno
tation/Trinotate/TransDecoder/true-unknown.fa.transdecoder-3.pep    

blastp -db ./localblastdb/localblastdb -query $PEP -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > WS-MS.outfmt
```     

This is awesome.... two really high confidence hits, and one less confident but still promising hit. 
1. STRG.12738.1.p1 as IPP2/D3Z3A0 - several peptide matches 
2. STRG.19052.1.p1 as PGAM - one of the high coverage hits    
3. STRG.8996.1.p1 as Takusan 

These are supported by PFAM and Uniref hits, have decently high transdecoder ORF scores also.      

In order to present this data, I will digest the primary peptide sequence on
Expasy PeptideMass service, and compare the digested fragments to the ones we
have recovered in our mass spec, to the ORFs predicted by Transdecoder. Also, I
will show the protein domains in SMART.     






