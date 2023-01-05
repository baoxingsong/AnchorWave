# MAF file to VCF file in the maize NAM

We introduce MAF file produced by AnchorWave converts to VCF file step by step. Maize NAM populations have 26 accessions and have abundant genetic diversity.

    software:
    AnchorWave: v1.0.1
    tassel: 5.2.82
    The Genome Analysis Toolkit (GATK): v4.2.6.1
    HTSJDK Version: 2.24.1
    Picard Version: 2.27.1
    Built for Spark Version: 2.4.5
    vt: v0.57721
    bcftools: 1.13 (using htslib 1.13+ds)
    igv: 2.15.4

## 1.1 Generate pair-wised genome alignments by AnchorWave

All the genome data were downloaded from maizegdb(<https://download.maizegdb.org/>). And CDS mapping were conducted via minimap2. The AnchorWave alignments were conducted using the following command. `anchorwave genoAli -i Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -as cds.fa -r B73.ref.fa -a B97.sam -ar B73.sam -s Zm-B97-REFERENCE-NAM-1.0.fa -n B97.anchors -o B97.maf -f B97.f.maf -w 38000 -fa3 200000 -B -6 -O1 -8 -E1 -2 -O2 -75 -E2 -1 -IV  >B97.log 2>&1`

## 1.2 MAF file convert to GVCF file

Convert MAF file to GVCF file by [tassel](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/UserInstructions/CreatePHG_step2_MAFToGVCFPluginDetails.md).

    #download tassel and reference genome.
    git clone https://bitbucket.org/tasseladmin/tassel-5-standalone.git 
    wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
    gunzip Zm-B73-REFERENCE-NAM-5.0.fa.gz

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/B97/B97.maf -sampleName B97_anchorwave -gvcfOutput /home/xuql/copyNAM/B97/B97ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/B97/B97_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/CML103/CML103.maf -sampleName CML103_anchorwave -gvcfOutput /home/xuql/copyNAM/CML103/CML103ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/CML103/CML103_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/CML228/CML228.maf -sampleName CML228_anchorwave -gvcfOutput /home/xuql/copyNAM/CML228/CML228ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/CML228/CML228_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/CML247/CML247.maf -sampleName CML247_anchorwave -gvcfOutput /home/xuql/copyNAM/CML247/CML247ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/CML247/CML247_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/CML277/CML277.maf -sampleName CML277_anchorwave -gvcfOutput /home/xuql/copyNAM/CML277/CML277ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/CML277/CML277_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/CML322/CML322.maf -sampleName CML322_anchorwave -gvcfOutput /home/xuql/copyNAM/CML322/CML322ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/CML322/CML322_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/CML333/CML333.maf -sampleName CML333_anchorwave -gvcfOutput /home/xuql/copyNAM/CML333/CML333ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/CML333/CML333_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/CML52/CML52.maf -sampleName CML52_anchorwave -gvcfOutput /home/xuql/copyNAM/CML52/CML52ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/CML52/CML52_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/CML69/CML69.maf -sampleName CML69_anchorwave -gvcfOutput /home/xuql/copyNAM/CML69/CML69ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/CML69/CML69_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/HP301/HP301.maf -sampleName HP301_anchorwave -gvcfOutput /home/xuql/copyNAM/HP301/HP301ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/HP301/HP301_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/Il14H/Il14H.maf -sampleName Il14H_anchorwave -gvcfOutput /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf -fillGaps false > /home/xuql/copyNAM/Il14H/Il14H_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/Ki11/Ki11.maf -sampleName Ki11_anchorwave -gvcfOutput /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/Ki11/Ki11_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/Ki3/Ki3.maf -sampleName Ki3_anchorwave -gvcfOutput /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/Ki3/Ki3_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/Ky21/Ky21.maf -sampleName Ky21_anchorwave -gvcfOutput /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/Ky21/Ky21_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/M162W/M162W.maf -sampleName M162W_anchorwave -gvcfOutput /home/xuql/copyNAM/M162W/M162WToB73.gvcf -fillGaps false > /home/xuql/copyNAM/M162W/M162W_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/M37W/M37W.maf -sampleName M37W_anchorwave -gvcfOutput /home/xuql/copyNAM/M37W/M37WToB73.gvcf -fillGaps false > /home/xuql/copyNAM/M37W/M37W_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/Mo18W/Mo18W.maf -sampleName Mo18W_anchorwave -gvcfOutput /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf -fillGaps false > /home/xuql/copyNAM/Mo18W/Mo18W_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/Ms71/Ms71.maf -sampleName Ms71_anchorwave -gvcfOutput /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/Ms71/Ms71_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/NC350/NC350.maf -sampleName NC350_anchorwave -gvcfOutput /home/xuql/copyNAM/NC350/NC350ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/NC350/NC350_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/NC358/NC358.maf -sampleName NC358_anchorwave -gvcfOutput /home/xuql/copyNAM/NC358/NC358ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/NC358/NC358_outputMafToGVCF.txt

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/Oh43/Oh43.maf -sampleName Oh43_anchorwave -gvcfOutput /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/Oh43/Oh43_outputMafToGVCF.txt 

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/Oh7B/Oh7B.maf -sampleName Oh7B_anchorwave -gvcfOutput /home/xuql/copyNAM/Oh7B/Oh7BToB73.gvcf -fillGaps false > /home/xuql/copyNAM/Oh7B/Oh7B_outputMafToGVCF.txt 

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/P39/P39.maf -sampleName P39_anchorwave -gvcfOutput /home/xuql/copyNAM/P39/P39ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/P39/P39_outputMafToGVCF.txt 

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/Tx303/Tx303.maf -sampleName Tx303_anchorwave -gvcfOutput /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/Tx303/Tx303_outputMafToGVCF.txt 

    /home/xuql/tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -mafFile /home/xuql/copyNAM/Tzi8/Tzi8.maf -sampleName Tzi8_anchorwave -gvcfOutput /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf -fillGaps false > /home/xuql/copyNAM/Tzi8/Tzi8_outputMafToGVCF.txt 

## 1.3 pre-treatment

GATK tool requires index and dictionary files for Reference Fasta files.

    sed -i 's/chr//g' Zm-B73-REFERENCE-NAM-5.0.fa 
    samtools faidx Zm-B73-REFERENCE-NAM-5.0.fa
    wget https://github.com/broadinstitute/picard/releases/download/2.26.10/picard.jar
    java -jar picard.jar CreateSequenceDictionary R=Zm-B73-REFERENCE-NAM-5.0.fa O=Zm-B73-REFERENCE-NAM-5.0.dict

The "GATK GenomicsDBImport" input file need be compressed and found index.

    bgzip -c  /home/xuql/copyNAM/B97/B97ToB73.gvcf > /home/xuql/copyNAM/B97/B97ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/CML103/CML103ToB73.gvcf > /home/xuql/copyNAM/CML103/CML103ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/CML228/CML228ToB73.gvcf > /home/xuql/copyNAM/CML228/CML228ToB73.gvcf.gz
    bgzip -c  /home/xuql/copyNAM/CML247/CML247ToB73.gvcf > /home/xuql/copyNAM/CML247/CML247ToB73.gvcf.gz
    bgzip -c  /home/xuql/copyNAM/CML277/CML277ToB73.gvcf > /home/xuql/copyNAM/CML277/CML277ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/CML322/CML322ToB73.gvcf > /home/xuql/copyNAM/CML322/CML322ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/CML333/CML333ToB73.gvcf > /home/xuql/copyNAM/CML333/CML333ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/CML52/CML52ToB73.gvcf > /home/xuql/copyNAM/CML52/CML52ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/CML69/CML69ToB73.gvcf > /home/xuql/copyNAM/CML69/CML69ToB73.gvcf.gz
    bgzip -c  /home/xuql/copyNAM/HP301/HP301ToB73.gvcf > /home/xuql/copyNAM/HP301/HP301ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf > /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf > /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf > /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf > /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/M162W/M162WToB73.gvcf > /home/xuql/copyNAM/M162W/M162WToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/M37W/M37WToB73.gvcf > /home/xuql/copyNAM/M37W/M37WToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf > /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf > /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/NC350/NC350ToB73.gvcf > /home/xuql/copyNAM/NC350/NC350ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/NC358/NC358ToB73.gvcf > /home/xuql/copyNAM/NC358/NC358ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf > /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/Oh7B/Oh7BToB73.gvcf > /home/xuql/copyNAM/Oh7B/Oh7BToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/P39/P39ToB73.gvcf > /home/xuql/copyNAM/P39/P39ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf > /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf.gz 
    bgzip -c  /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf > /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf.gz

    tabix -p vcf /home/xuql/copyNAM/B97/B97ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/CML247/CML247ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/CML333/CML333ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/HP301/HP301ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/M37W/M37WToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/NC350/NC350ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/Oh7B/Oh7BToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/CML103/CML103ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/CML277/CML277ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/CML52/CML52ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/NC358/NC358ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/P39/P39ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/CML228/CML228ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/CML322/CML322ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/CML69/CML69ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/M162W/M162WToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf.gz 

## 1.4 Merge GVCF files into a GATK database.

We conducted this analysis from chr 1 to chr 10 separately. Please note that the longest indel could not over 9101264 we tested, else,the next step would report error. We have posed the question to [GATK#7976](https://github.com/broadinstitute/gatk/issues/7976). To check indel length distribution, you could use "gatk LeftAlignAndTrimVariants". We filtered \>9101264 indels that only exist chromosome 9 and 10.

check out and filter super-indels procedure:

    #for indel how to sort by python and the part results is shown.
    with open('./job.407755.err') as f:
        lines = f.readlines()
    new_lines=[]
    for line in lines:
        if "Set" in line:
            new_lines.append(line)
            
    sorted(new_lines, 
           key=lambda x: int(x.strip().split('=')[-1]),
           reverse=True)

    #thoese content output from all chromosome by "gatk LeftAlignAndTrimVariants" and then sorted by python. at last,check the initial output content to identify the super-indel belong to individuals.
    ['10:56:14.335 INFO  LeftAlignAndTrimVariants - Indel is too long (34461688) at position 9:3695105; skipping that record. Set --max-indel-length >= 34461688\n',
     '10:56:30.429 INFO  LeftAlignAndTrimVariants - Indel is too long (10668738) at position 10:33212598; skipping that record. Set --max-indel-length >= 10668738\n',
     '10:56:28.937 INFO  LeftAlignAndTrimVariants - Indel is too long (9101264) at position 10:14179; skipping that record. Set --max-indel-length >= 9101264\n',
     '10:56:30.038 INFO  LeftAlignAndTrimVariants - Indel is too long (7918835) at position 10:22996027; skipping that record. Set --max-indel-length >= 7918835\n',
     '11:31:49.968 INFO  LeftAlignAndTrimVariants - Indel is too long (7154442) at position 6:16715313; skipping that record. Set --max-indel-length >= 7154442\n',

    #acquire the line number of super-indel and delete it.
    less Oh7BToB73.gvcf | awk '$2=="3695105"{printf("%5d\t%s\n", NR, $0)}' >aaa.txt 
    less Oh7BToB73.gvcf | sed '32542168d' > aOh7BToB73.gvcf 
    # compress and index for "gatk GenomicsDBImport"
    bgzip -c  /home/xuql/copyNAM/Oh7B/aOh7BToB73.gvcf > /home/xuql/copyNAM/Oh7B/aOh7BToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/Oh7B/aOh7BToB73.gvcf.gz

    less Oh7BToB73.gvcf | awk '$2=="33212598"{printf("%5d\t%s\n", NR, $0)}' >bbb.txt
    less Oh7BToB73.gvcf | sed '35176648d' > bOh7BToB73.gvcf
    bgzip -c  /home/xuql/copyNAM/Oh7B/bOh7BToB73.gvcf > /home/xuql/copyNAM/Oh7B/bOh7BToB73.gvcf.gz 
    tabix -p vcf /home/xuql/copyNAM/Oh7B/bOh7BToB73.gvcf.gz

    gatk --java-options "-Xmx256g -Xms256g" GenomicsDBImport \
                    -V /home/xuql/copyNAM/B97/B97ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML247/CML247ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML333/CML333ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/HP301/HP301ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M37W/M37WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC350/NC350ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh7B/Oh7BToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML103/CML103ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML277/CML277ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML52/CML52ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC358/NC358ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/P39/P39ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML228/CML228ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML322/CML322ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML69/CML69ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M162W/M162WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf.gz \
            --batch-size 5 \
          --genomicsdb-workspace-path /home/xuql/copyNAM/NAM_out_gatk1 \
          --genomicsdb-segment-size 1048576000 --genomicsdb-vcf-buffer-size 10000000000 -L 1
          #Elapsed time: 105.02 minutes. Runtime.totalMemory()=274877906944

    gatk --java-options "-Xmx256g -Xms256g" GenomicsDBImport \
                    -V /home/xuql/copyNAM/B97/B97ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML247/CML247ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML333/CML333ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/HP301/HP301ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M37W/M37WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC350/NC350ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh7B/Oh7BToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML103/CML103ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML277/CML277ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML52/CML52ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC358/NC358ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/P39/P39ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML228/CML228ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML322/CML322ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML69/CML69ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M162W/M162WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf.gz \
            --batch-size 1 \
          --genomicsdb-workspace-path /home/xuql/copyNAM/NAM_out_gatk2 \
          --genomicsdb-segment-size 10485760 --genomicsdb-vcf-buffer-size 100000000 -L 2
            #Elapsed time: 76.68 minutes.  Runtime.totalMemory()=274877906944

    gatk --java-options "-Xmx256g -Xms256g" GenomicsDBImport \
                    -V /home/xuql/copyNAM/B97/B97ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML247/CML247ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML333/CML333ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/HP301/HP301ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M37W/M37WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC350/NC350ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh7B/Oh7BToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML103/CML103ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML277/CML277ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML52/CML52ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC358/NC358ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/P39/P39ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML228/CML228ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML322/CML322ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML69/CML69ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M162W/M162WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf.gz \
            --batch-size 1 \
          --genomicsdb-workspace-path /home/xuql/copyNAM/NAM_out_gatk3 \
          --genomicsdb-segment-size 10485760 --genomicsdb-vcf-buffer-size 100000000 -L 3
            #Elapsed time: 70.12 minutes. Runtime.totalMemory()=274877906944

    gatk --java-options "-Xmx256g -Xms256g" GenomicsDBImport \
                    -V /home/xuql/copyNAM/B97/B97ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML247/CML247ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML333/CML333ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/HP301/HP301ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M37W/M37WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC350/NC350ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh7B/Oh7BToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML103/CML103ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML277/CML277ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML52/CML52ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC358/NC358ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/P39/P39ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML228/CML228ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML322/CML322ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML69/CML69ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M162W/M162WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf.gz \
            --batch-size 1 \
          --genomicsdb-workspace-path /home/xuql/copyNAM/NAM_out_gatk4 \
          --genomicsdb-segment-size 10485760 --genomicsdb-vcf-buffer-size 100000000 -L 4
            #Elapsed time: 68.47 minutes. Runtime.totalMemory()=274877906944

    gatk --java-options "-Xmx256g -Xms256g" GenomicsDBImport \
                    -V /home/xuql/copyNAM/B97/B97ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML247/CML247ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML333/CML333ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/HP301/HP301ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M37W/M37WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC350/NC350ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh7B/Oh7BToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML103/CML103ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML277/CML277ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML52/CML52ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC358/NC358ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/P39/P39ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML228/CML228ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML322/CML322ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML69/CML69ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M162W/M162WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf.gz \
            --batch-size 1 \
          --genomicsdb-workspace-path /home/xuql/copyNAM/NAM_out_gatk5 \
          --genomicsdb-segment-size 10485760 --genomicsdb-vcf-buffer-size 100000000 -L 5
            # Elapsed time: 64.29 minutes. Runtime.totalMemory()=274877906944

    gatk --java-options "-Xmx256g -Xms256g" GenomicsDBImport \
                    -V /home/xuql/copyNAM/B97/B97ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML247/CML247ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML333/CML333ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/HP301/HP301ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M37W/M37WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC350/NC350ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh7B/Oh7BToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML103/CML103ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML277/CML277ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML52/CML52ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC358/NC358ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/P39/P39ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML228/CML228ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML322/CML322ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML69/CML69ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M162W/M162WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf.gz \
            --batch-size 1 \
          --genomicsdb-workspace-path /home/xuql/copyNAM/NAM_out_gatk6 \
          --genomicsdb-segment-size 10485760 --genomicsdb-vcf-buffer-size 100000000 -L 6
            #Elapsed time: 53.32 minutes. Runtime.totalMemory()=274877906944

    gatk --java-options "-Xmx256g -Xms256g" GenomicsDBImport \
                    -V /home/xuql/copyNAM/B97/B97ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML247/CML247ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML333/CML333ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/HP301/HP301ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M37W/M37WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC350/NC350ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh7B/Oh7BToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML103/CML103ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML277/CML277ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML52/CML52ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC358/NC358ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/P39/P39ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML228/CML228ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML322/CML322ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML69/CML69ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M162W/M162WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf.gz \
            --batch-size 1 \
          --genomicsdb-workspace-path /home/xuql/copyNAM/NAM_out_gatk7 \
          --genomicsdb-segment-size 10485760 --genomicsdb-vcf-buffer-size 100000000 -L 7
            #Elapsed time: 57.26 minutes.  Runtime.totalMemory()=274877906944

    gatk --java-options "-Xmx256g -Xms256g" GenomicsDBImport \
                    -V /home/xuql/copyNAM/B97/B97ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML247/CML247ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML333/CML333ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/HP301/HP301ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M37W/M37WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC350/NC350ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh7B/Oh7BToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML103/CML103ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML277/CML277ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML52/CML52ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC358/NC358ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/P39/P39ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML228/CML228ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML322/CML322ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML69/CML69ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M162W/M162WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf.gz \
            --batch-size 1 \
          --genomicsdb-workspace-path /home/xuql/copyNAM/NAM_out_gatk8 \
          --genomicsdb-segment-size 10485760 --genomicsdb-vcf-buffer-size 100000000 -L 8
          #Elapsed time: 52.72 minutes. Runtime.totalMemory()=274877906944

    gatk --java-options "-Xmx128g -Xms5g" GenomicsDBImport \
                    -V /home/xuql/copyNAM/B97/B97ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML247/CML247ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML333/CML333ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/HP301/HP301ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M37W/M37WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC350/NC350ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh7B/aOh7BToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML103/CML103ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML277/CML277ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML52/CML52ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC358/NC358ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/P39/P39ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML228/CML228ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML322/CML322ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML69/CML69ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M162W/M162WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf.gz \
            --batch-size 1 \
          --genomicsdb-workspace-path /home/xuql/copyNAM/NAM_out_gatk9 \
          --genomicsdb-segment-size 10485760 --genomicsdb-vcf-buffer-size 100000000 -L 9
          #Elapsed time: 54.09 minutes.Runtime.totalMemory()=5368709120

    gatk --java-options "-Xmx128g -Xms5g" GenomicsDBImport \
                    -V /home/xuql/copyNAM/B97/B97ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML247/CML247ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML333/CML333ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/HP301/HP301ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki3/Ki3ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M37W/M37WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC350/NC350ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh7B/bOh7BToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tzi8/Tzi8ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML103/CML103ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML277/CML277ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML52/CML52ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Il14H/Il14HToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ky21/Ky21ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Mo18W/Mo18WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/NC358/NC358ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/P39/P39ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML228/CML228ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML322/CML322ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/CML69/CML69ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ki11/Ki11ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/M162W/M162WToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Ms71/Ms71ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Oh43/Oh43ToB73.gvcf.gz \
                    -V /home/xuql/copyNAM/Tx303/Tx303ToB73.gvcf.gz \
            --batch-size 1 \
          --genomicsdb-workspace-path /home/xuql/copyNAM/NAM_out_gatk10 \
          --genomicsdb-segment-size 10485760 --genomicsdb-vcf-buffer-size 100000000 -L 10  
            #Elapsed time: 44.05 minutes.Runtime.totalMemory()=5368709120

## 1.5 Vairant calling

Please note that the max indel don't over 9101264 we tested, for test max indel length you can do "gatk LeftAlignAndTrimVariants" at default parameters for all chromosome and sorted it by python. and we have posed the question to GATK[#7976](https://github.com/broadinstitute/gatk/issues/7976).

    gatk --java-options "-Xmx50g" GenotypeGVCFs -R /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -stand-call-conf 0 -ploidy 1 -V gendb:///home/xuql/copyNAM/NAM_out_gatk1 -O /home/xuql/copyNAM/gatk1.vcf.gz --cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000
    #Processed 251973478 total variants in 581.3 minutes. Elapsed time: 581.36 minutes. Runtime.totalMemory()=6694109184

    gatk --java-options "-Xmx50g" GenotypeGVCFs -R /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -stand-call-conf 0 -ploidy 1 -V gendb:///home/xuql/copyNAM/NAM_out_gatk2 -O /home/xuql/copyNAM/gatk2.vcf.gz --cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000
    # Processed 202340450 total variants in 471.8 minutes.Elapsed time: 472.05 minutes.Runtime.totalMemory()=2147483648

    gatk --java-options "-Xmx50g" GenotypeGVCFs -R /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -stand-call-conf 0 -ploidy 1 -V gendb:///home/xuql/copyNAM/NAM_out_gatk3 -O /home/xuql/copyNAM/gatk3.vcf.gz --cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000
    # Processed 188148248 total variants in 444.9 minutes.Elapsed time: 445.13 minutes.Runtime.totalMemory()=2147483648

    gatk --java-options "-Xmx50g" GenotypeGVCFs -R /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -stand-call-conf 0 -ploidy 1 -V gendb:///home/xuql/copyNAM/NAM_out_gatk4 -O /home/xuql/copyNAM/gatk4.vcf.gz --cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000
    #Processed 188872000 total variants in 432.5 minutes. Elapsed time: 432.67 minutes.Runtime.totalMemory()=2147483648

    gatk --java-options "-Xmx50g" GenotypeGVCFs -R /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -stand-call-conf 0 -ploidy 1 -V gendb:///home/xuql/copyNAM/NAM_out_gatk5 -O /home/xuql/copyNAM/gatk5.vcf.gz --cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000
    #Elapsed time: 434.09 minutes.Runtime.totalMemory()=2147483648
     
    gatk --java-options "-Xmx50g" GenotypeGVCFs -R /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -stand-call-conf 0 -ploidy 1 -V gendb:///home/xuql/copyNAM/NAM_out_gatk6 -O /home/xuql/copyNAM/gatk6.vcf.gz --cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000
    # Processed 148426585 total variants in 338.4 minutes.Elapsed time: 338.57 minutes. Runtime.totalMemory()=7088373760

    gatk --java-options "-Xmx50g" GenotypeGVCFs -R /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -stand-call-conf 0 -ploidy 1 -V gendb:///home/xuql/copyNAM/NAM_out_gatk7 -O /home/xuql/copyNAM/gatk7.vcf.gz --cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000
    #Processed 148023151 total variants in 369.9 minutes.Elapsed time: 370.07 minutes. Runtime.totalMemory()=2147483648

    gatk --java-options "-Xmx50g" GenotypeGVCFs -R /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -stand-call-conf 0 -ploidy 1 -V gendb:///home/xuql/copyNAM/NAM_out_gatk8 -O /home/xuql/copyNAM/gatk8.vcf.gz --cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000
    #Processed 147858803 total variants in 355.4 minutes.Elapsed time: 355.56 minutes.Runtime.totalMemory()=7214202880

    #to make sure "GenotypeGVCFs" can implement we delete the indels more than 10M that only included in cheomosome 9 and chromosome 10.

    gatk --java-options "-Xmx100g" GenotypeGVCFs -R /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -stand-call-conf 0 -ploidy 1 -V gendb:///home/xuql/copyNAM/NAM_out_gatk9 -O /home/xuql/copyNAM/gatk9.vcf.gz --cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000
    #Elapsed time: 320.96 minutes.Runtime.totalMemory()=3724541952

    gatk --java-options "-Xmx200g" GenotypeGVCFs -R /home/xuql/copyNAM/B73/Zm-B73-REFERENCE-NAM-5.0.fa -stand-call-conf 0 -ploidy 1 -V gendb:///home/xuql/copyNAM/NAM_out_gatk10 -O /home/xuql/copyNAM/gatk10.vcf.gz --cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000
    #Elapsed time: 291.33 minutes.Runtime.totalMemory()=7079985152

## 1.6 Variation normalization

    # We used vt to noamalize the indel and bcftools to merge duplication at the same varition loci.

    /home/ywt/vt/vt normalize gatk1.vcf.gz -r /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -o gatk1bcfvt.vcf.gz 
    bcftools norm  -f /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -m +both gatk1bcfvt.vcf.gz -Oz -o gatk1bcftools.vcf.gz 
    bcftools index -t gatk1bcfvt.vcf.gz
    bcftools index -t gatk1bcftools.vcf.gz

    /home/ywt/vt/vt normalize gatk2.vcf.gz -r /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -o gatk2bcfvt.vcf.gz 
    bcftools norm  -f /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -m +both gatk2bcfvt.vcf.gz -Oz -o gatk2bcftools.vcf.gz 
    bcftools index -t gatk2bcfvt.vcf.gz
    bcftools index -t gatk2bcftools.vcf.gz

    /home/ywt/vt/vt normalize gatk3.vcf.gz -r /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -o gatk3bcfvt.vcf.gz 
    bcftools norm  -f /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -m +both gatk3bcfvt.vcf.gz -Oz -o gatk3bcftools.vcf.gz 
    bcftools index -t gatk3bcfvt.vcf.gz
    bcftools index -t gatk3bcftools.vcf.gz

    /home/ywt/vt/vt normalize gatk4.vcf.gz -r /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -o gatk4bcfvt.vcf.gz 
    bcftools norm  -f /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -m +both gatk4bcfvt.vcf.gz -Oz -o gatk4bcftools.vcf.gz 
    bcftools index -t gatk4bcfvt.vcf.gz
    bcftools index -t gatk4bcftools.vcf.gz

    /home/ywt/vt/vt normalize gatk5.vcf.gz -r /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -o gatk5bcfvt.vcf.gz 
    bcftools norm  -f /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -m +both gatk5bcfvt.vcf.gz -Oz -o gatk5bcftools.vcf.gz 
    bcftools index -t gatk5bcfvt.vcf.gz
    bcftools index -t gatk5bcftools.vcf.gz

    /home/ywt/vt/vt normalize gatk6.vcf.gz -r /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -o gatk6bcfvt.vcf.gz 
    bcftools norm  -f /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -m +both gatk6bcfvt.vcf.gz -Oz -o gatk6bcftools.vcf.gz 
    bcftools index -t gatk6bcfvt.vcf.gz
    bcftools index -t gatk6bcftools.vcf.gz

    /home/ywt/vt/vt normalize gatk7.vcf.gz -r /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -o gatk7bcfvt.vcf.gz 
    bcftools norm  -f /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -m +both gatk7bcfvt.vcf.gz -Oz -o gatk7bcftools.vcf.gz 
    bcftools index -t gatk7bcfvt.vcf.gz
    bcftools index -t gatk7bcftools.vcf.gz

    bcftools index -t gatk8.vcf.gz
    /home/ywt/vt/vt normalize gatk8.vcf.gz -r /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -o gatk8bcfvt.vcf.gz 
    bcftools norm -f /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -m +both gatk8bcfvt.vcf.gz -Oz -o gatk8bcftools.vcf.gz
    bcftools index -t gatk8bcfvt.vcf.gz
    bcftools index -t gatk8bcftools.vcf.gz

    /home/ywt/vt/vt normalize gatk9.vcf.gz -r /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -o gatk9bcfvt.vcf.gz 
    bcftools norm -f /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -m +both gatk9bcfvt.vcf.gz -Oz -o gatk9bcftools.vcf.gz
    bcftools index -t gatk9bcfvt.vcf.gz
    bcftools index -t gatk9bcftools.vcf.gz

    /home/ywt/vt/vt normalize gatk10.vcf.gz -r /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -o gatk10bcfvt.vcf.gz 
    bcftools norm -f /media/ywt/14T1/NAManchorwave/B73/Zm-B73-REFERENCE-NAM-5.0.fa -m +both gatk10bcfvt.vcf.gz -Oz -o gatk10bcftools.vcf.gz
    bcftools index -t gatk10bcfvt.vcf.gz
    bcftools index -t gatk10bcftools.vcf.gz

## 1.7 IGV visualize
The genome alignments could be converted to cram/bam files

```maf-convert sam B97.maf > B97.sam
sed -r -i 's/[0-9]+H//g' B97.sam &
samtools view -O CRAM --threads 80 --reference Zm-B73-REFERENCE-NAM-5.0.fa B97.sam | samtools sort --threads 30 -O CRAM - > B97.cram
samtools index B97.cram
```
cram files and the vcf file could be loaded into IGV to check the alginments and variant calling manually.

If you have any questions, please feel free to contact Shuai Wang E-mail: [shuaiwang2\@163.com](mailto:shuaiwang2@163.com){.email} or wechat:shuaiwang2022.
