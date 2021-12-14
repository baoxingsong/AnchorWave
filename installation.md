# Installation
AnchorWave uses the SIMD instructions to speed up the dynamic programming algorithm. Specific functions have been implemented for SSE2, SSE4.1, AVX2 and AVX512 instruction sets.  
Empirically, the time cost: SSE2 > SSE4.1 > AVX2 > AVX512.  
The conda version compile all of them and pick the correct version to run.  
To check what CPU instructions are supported by your machine, you could run this command:
```cat /proc/cpuinfo | grep "flags" | uniq```  
By default, we assume the machine supports SSE4.1.  

## If you are using machine with AVX512 and would like to take the advantage of that
Clone the repository, and replace the default CMakeLists.txt with the avx512 one.
```
git clone https://github.com/baoxingsong/anchorwave.git
cd anchorwave
mv CMakeLists_avx512.txt CMakeLists.txt
cmake ./
make
```


## If you are using machine with AVX2 and would like to take the advantage of AVX2
Clone the repository, and replace the default CMakeLists.txt with the avx2 one.
```
git clone https://github.com/baoxingsong/anchorwave.git
cd anchorwave
mv CMakeLists_avx2.txt CMakeLists.txt
cmake ./
make
```


## If you are using old machine with SSE2 but not SSE4.1

Clone the repository, and replace the default CMakeLists.txt with the sse2 one.
```
git clone https://github.com/baoxingsong/anchorwave.git
cd anchorwave
mv CMakeLists_sse2.txt CMakeLists.txt
cmake ./
make
```

## If you are using very old CPU or CPU with other architecture (i.e. ARM)
AnchorWave was not tested on that kind of platforms. High likely would not work properly.

# Time cost comparison using different installation ways
## Hardware
We compiled AnchorWave in 4 different ways on a machine that having two Intel(R) Xeon(R) Gold 6238 CPUs and 512Gb RAM.
## Data and method
We aligned the maize Mo17 CAU assembly against the B73 V4 reference assembly. Input data were prepared using the following commands:
```
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-34/gff3/zea_mays/Zea_mays.AGPv4.34.gff3.gz
gunzip Zea_mays.AGPv4.34.gff3.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-34/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
gunzip Zea_mays.AGPv4.dna.toplevel.fa.gz
wget https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-1.0/Zm-Mo17-REFERENCE-CAU-1.0.fa.gz
gunzip Zm-Mo17-REFERENCE-CAU-1.0.fa.gz
sed -i 's/>chr/>/g' Zm-Mo17-REFERENCE-CAU-1.0.fa

./anchorwave_avx512 gff2seq -r Zea_mays.AGPv4.dna.toplevel.fa -i Zea_mays.AGPv4.34.gff3 -o cds.fa

minimap2 -x splice -t 80 -k 12 -a -p 0.4 -N 20 Zm-Mo17-REFERENCE-CAU-1.0.fa cds.fa > cds.sam
minimap2 -x splice -t 80 -k 12 -a -p 0.4 -N 20 Zea_mays.AGPv4.dna.toplevel.fa cds.fa > ref.sam
```
Then we run the following 4 commands on the same computer simultaneously:
```
/usr/bin/time ./anchorwave_sse2 genoAli -i Zea_mays.AGPv4.34.gff3 -as cds.fa -r Zea_mays.AGPv4.dna.toplevel.fa -a cds.sam -ar ref.sam -s Zm-Mo17-REFERENCE-CAU-1.0.fa -n sse2.anchors -o sse2.anchorwave.maf -f sse2.anchorwave.f.maf -w 38000 -fa3 200000 -B -6 -O1 -8 -E1 -2 -O2 -75 -E2 -1 -t 2 -IV >sse2.log 2>&1
/usr/bin/time ./anchorwave_sse4.1 genoAli -i Zea_mays.AGPv4.34.gff3 -as cds.fa -r Zea_mays.AGPv4.dna.toplevel.fa -a cds.sam -ar ref.sam -s Zm-Mo17-REFERENCE-CAU-1.0.fa -n sse4.1.anchors -o sse4.1.anchorwave.maf -f sse4.1.anchorwave.f.maf -w 38000 -fa3 200000 -B -6 -O1 -8 -E1 -2 -O2 -75 -E2 -1 -t 2 -IV >sse4.1.log 2>&1
/usr/bin/time ./anchorwave_avx2 genoAli -i Zea_mays.AGPv4.34.gff3 -as cds.fa -r Zea_mays.AGPv4.dna.toplevel.fa -a cds.sam -ar ref.sam -s Zm-Mo17-REFERENCE-CAU-1.0.fa -n avx2.anchors -o avx2.anchorwave.maf -f avx2.anchorwave.f.maf -w 38000 -fa3 200000 -B -6 -O1 -8 -E1 -2 -O2 -75 -E2 -1 -t 2 -IV >avx2.log 2>&1
/usr/bin/time ./anchorwave_avx512 genoAli -i Zea_mays.AGPv4.34.gff3 -as cds.fa -r Zea_mays.AGPv4.dna.toplevel.fa -a cds.sam -ar ref.sam -s Zm-Mo17-REFERENCE-CAU-1.0.fa -n avx512.anchors -o avx512.anchorwave.maf -f avx512.anchorwave.f.maf -w 38000 -fa3 200000 -B -6 -O1 -8 -E1 -2 -O2 -75 -E2 -1 -t 2 -IV >avx512.log 2>&1
```
Compare the output alignments:
```
maf-convert sam sse2.anchorwave.maf > sse2.anchorwave.sam
maf-convert sam sse4.1.anchorwave.maf > sse4.1.anchorwave.sam
maf-convert sam avx2.anchorwave.maf > avx2.anchorwave.sam
maf-convert sam avx512.anchorwave.maf > avx512.anchorwave.sam

sed -i -r 's/[0-9]+H//g' sse2.anchorwave.sam
sed -i -r 's/[0-9]+H//g' sse4.1.anchorwave.sam
sed -i -r 's/[0-9]+H//g' avx2.anchorwave.sam
sed -i -r 's/[0-9]+H//g' avx512.anchorwave.sam

samtools view -O CRAM --reference Zea_mays.AGPv4.dna.toplevel.fa sse2.anchorwave.sam | samtools sort -O CRAM - > sse2.anchorwave.cram
samtools view -O CRAM --reference Zea_mays.AGPv4.dna.toplevel.fa sse4.1.anchorwave.sam | samtools sort -O CRAM - > sse4.1.anchorwave.cram
samtools view -O CRAM --reference Zea_mays.AGPv4.dna.toplevel.fa avx2.anchorwave.sam | samtools sort -O CRAM - > avx2.anchorwave.cram
samtools view -O CRAM --reference Zea_mays.AGPv4.dna.toplevel.fa avx512.anchorwave.sam | samtools sort -O CRAM - > avx512.anchorwave.cram

samtools view -O SAM -o sse2.anchorwave.sorted.sam sse2.anchorwave.cram
samtools view -O SAM -o sse4.1.anchorwave.sorted.sam sse4.1.anchorwave.cram
samtools view -O SAM -o avx2.anchorwave.sorted.sam avx2.anchorwave.cram
samtools view -O SAM -o avx512.anchorwave.sorted.sam avx512.anchorwave.cram

diff sse4.1.anchorwave.sorted.sam avx2.anchorwave.sorted.sam
diff sse4.1.anchorwave.sorted.sam avx512.anchorwave.sorted.sam
diff sse2.anchorwave.sorted.sam avx512.anchorwave.sorted.sam
```
Different compiling ways generated identical alignments
## RAM and wall time cost comparison
|             | peak memory (Gb) | wall time   |
| ----------- | ----------- | ----------- | 
| SSE2        | 134.9        | 38:20:02 |
| SSE4.1      | 132.5        | 35:52:18 |
| AVX2        | 142.4        | 35:04:26 |
| AVX512      | 133.2        | 32:01:02 |

*the cost comparison was performed for once (did not perform multiple times and calculate the average)
