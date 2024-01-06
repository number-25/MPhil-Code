# Document outline     
Herein, I will explore the thesis that long-reads offer advantages in pseudogene
discovery due to their extended coverage (read length) of discriminatory
regions in pseudogenes. As pseudogenes ordinarily have high sequence similarity
to their parent genes, long reads may allow us to cross unique regions in the
pseudogene with single, unfragmented reads, boosting information content, thus
allowing for effective discrimination between pseudogenes and their bona-fide
gene.      



* If the sequence is differentiated by only a single or a few point mutations/changes, if these are manually made to the sequence, does the fastq map to the parent gene instead and not the pseudogene?    


Example codes for doing this     
```
samtools view RC-EXT.bam chr10:108045776-108046446 -o tmp.bam

samtools fastq tmp.bam | less -S

samtools fasta tmp.bam | less -S
```    

Create a .bam file from the positive intersects of the dRNA pseudogene set with the in-house short-read data.     

```
less -S *.bed | cut -f1,2,3,4,5,6,7,8,9,10 | samtools view -@ 10 -L -
/media/labpc/Disk-3/Bioinformatics/tmp-analysis-output-data/Direct-RNA-July2019-Garvan/Analysis/pending-upload/HISAT2/merged/merged/out.sorted.bam
> PC-NC_PS_Short-read-Intersect.bam     
```     

Generate some stats for this newly generated .bam file.    
`NanoStat -t8 --bam *.bam > nano-bam.stats`     

Create a few bam plots also.     
`NanoPlot -t 12 --color yellow --bam PC-NC_PS_Short-read-Intersect.bam -o bamplots`       

> General summary:        
Average percent identity:        99.8
Mean read length:              157.6
Mean read quality:              38.9
Median percent identity:       100.0
Median read length:            150.5
Median read quality:            40.5
Number of reads:             2,294.0
Read length N50:               168.0
Total bases:               361,582.0

* High read quality on average - 150bp read lengths, pretty nice.     


Do as above, but with the pseudogene long-read data .bam - this will give an idea of the quality scores of the fastq reads.    

Answering the question: Can novel annotations be reliably made even amongst the known indel-error rate that nanopore carries? The fact that we have
uncovered a few should be a strong affirmative.

`NanoStat -t 8 --bam NC-PC-Pseudogene-set-transcript.bam > nano-bam_dRNA.stats`      

> General summary:        
Average percent identity:         88.9
Mean read length:               464.7
Mean read quality:                8.6
Median percent identity:         89.8
Median read length:             370.0
Median read quality:              8.6
Number of reads:              3,666.0
Read length N50:                486.0
Total bases:              1,703,748.0
Total bases aligned:      1,547,571.0
Number, percentage and megabases of reads above quality cutoffs
>Q5:    3610 (98.5%) 1.7Mb
>Q7:    3070 (83.7%) 1.5Mb
>Q10:   753 (20.5%) 0.5Mb
>Q12:   56 (1.5%) 0.0Mb
>Q15:   2 (0.1%) 0.0Mb
Top 5 highest mean basecall quality scores and their read lengths
1:      15.2 (149; aebde6ee-3ec2-4d3e-b110-ac9596e944b5)
2:      15.1 (175; 9c18af10-0221-42a0-af81-c292079a0f39)
3:      14.3 (108; fbcd9b91-a465-43f6-aa7e-66582cdaff56)
4:      13.7 (793; 4cc00163-2877-4ef6-853c-3d38553f4a20)
5:      13.5 (283; 554474d2-0d84-4104-91e1-d152344db1b2)
Top 5 longest reads and their mean basecall quality score
1:      5447 (11.1; 7b31c91a-43a0-4318-a17b-7ca3f6ae6b05)
2:      5068 (8.7; 00bbfea2-0569-4c46-9c0b-b902fbcd9395)
3:      4176 (5.1; 61dd0625-5378-4cf3-90b4-ba39b97665eb)
4:      4099 (10.3; 44f91e1f-cfb1-4755-a47f-2a652e98c23b)
5:      3881 (11.9; 9d19c5af-51bf-4831-a841-e6800e3a6370)

* What's interesting here is: 1) the highest quality scores have read lengths
virtually in-line with short-read technologies, although at much lower
quality scores. 2) The longest reads are magnitudes longer than that offered
by short-read tech. 3) Overall, average read lengths are quite, and quality
scores and quite comparatively low, yet novel pseudogene discovery and
annotation is still possible.       

Generate some .bam plots.    
`NanoPlot -t 12 --color green --bam NC-PC-Pseudogene-set-transcript.bam -o bamplots`      

