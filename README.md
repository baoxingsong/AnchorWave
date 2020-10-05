# progli &middot; [![license][license-badge]][license]

Using the alignment of conversed sequence (currently, only CDS being used), progli try to detect conserved collinear blocks between two genomes. Within each collinear block, progli use conserved sequence as anchors and cut long sequence into short fragments and perform global sequence alignment for each block.
If the block length is longer than a preset window size, the output alignment maybe problematic. While a large window size would cost more memories.



## Install
### Dependencies
GNU GCC >=6.0 \
Cmake >= 3.0 \
minimap2
```
git clone https://github.com/baoxingsong/proali.git
cd proali
cmake ./
make
```

this command will generate an executable file named proali

### proportional genome alignment 
data: maize B73 genome sequence and GFF3 annotation file from http://plants.ensembl.org/Zea_mays/Info/Index
and sorghum genome sequence from https://plants.ensembl.org/Sorghum_bicolor/Info/Index

extract CDS of reference specie and map the to the query genome using minimap2 
```
proali gff2seq -i Zea_mays.AGPv4.34.gff3 -r Zea_mays.AGPv4.dna.toplevel.fa  -o cds.fa
minimap2 -ax splice -t 90 -a -uf -p 0.4 -C5 Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa ../primary_cds.fasta >cds.sam
```

```
proali proali -i Zea_mays.AGPv4.34.gff3 -r Zea_mays.AGPv4.dna.toplevel.fa -a cds.sam -s Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa -w 20000 -R 1 -Q 2 -O 0 -D 32 -o alignment -f
```

### Variant calling for different accessions
data: Arabidopsis thaliana Col-o reference genome and GFF3 annotation file from https://www.arabidopsis.org/
Arabidopsis thaliana Ler-0 accession assembly from http://www.pnas.org/content/113/28/E4052
 
```
proali gff2seq -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -o cds.fa
minimap2 -ax splice -t 10 -a -uf -p 0.4 -C5 ler.fa cds.fa >cds.sam
```
please make sure the chromosomes from reference genome and query genomes were named in the same way. Chromosomes with the same names would be aligned
```
grep ">" ler.fa
grep ">" Col.fa
```
perform alignmetn and variant calling
```
proali genoAli -i TAIR10_GFF3_genes_no_UM.gff -r tair10.fa -a cds.sam -s ler.fa  -f fragmentedAlignment.maf -v variant.vcf -m alignment -o anchorBlocks.tsv
```

### Contact
Bug report? Any question? Any suggestion? Any requirement? Want to cooperate?\
Please feel free to send E-mail to songbaoxing168@163.com

### Founding
This work is funded by NSF #1822330

### Citation
The manuscript is under preparation

[license]: ./LICENSE
[license-badge]: https://img.shields.io/badge/license-MIT-blue.svg
