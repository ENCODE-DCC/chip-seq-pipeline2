
WORKDIR=/mnt/data/pipeline_test_samples/chipseq/ENCSR936XTK/fastq_chr21
SAMPLE=ENCSR936XTK
CHR=chr21
BAM1=/srv/scratch/shared/surya/leepc12/run/NatProtocol_TestCase/ENCSR936XTK/PE/align/rep1/rep1-R1.PE2SE.bam
BAM2=/srv/scratch/shared/surya/leepc12/run/NatProtocol_TestCase/ENCSR936XTK/PE/align/rep2/rep2-R1.PE2SE.bam
SMALL_BAM1=$SAMPLE.$CHR.rep1.bam
SMALL_BAM2=$SAMPLE.$CHR.rep2.bam
SMALL_NMSRT_BAM1=$SAMPLE.$CHR.rep1.nmsrt.bam
SMALL_NMSRT_BAM2=$SAMPLE.$CHR.rep2.nmsrt.bam
SMALL_FQ1_R1=$SAMPLE.$CHR.rep1-R1.fastq
SMALL_FQ1_R2=$SAMPLE.$CHR.rep1-R2.fastq
SMALL_FQ2_R1=$SAMPLE.$CHR.rep2-R1.fastq
SMALL_FQ2_R2=$SAMPLE.$CHR.rep2-R2.fastq
mkdir -p $WORKDIR && cd $WORKDIR
samtools view -b $BAM1 "$CHR" > $SMALL_BAM1
samtools view -b $BAM2 "$CHR" > $SMALL_BAM2
samtools sort -n $SMALL_BAM1 $SMALL_NMSRT_BAM1
samtools sort -n $SMALL_BAM2 $SMALL_NMSRT_BAM2
#samtools bamshuf -uO $SMALL_NMSRT_BAM1.bam tmp-prefix | samtools bam2fq -s se.fq.gz - | bwa mem -p /mnt/data/pipeline_genome_data/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta - > tmp.fq1
#samtools bamshuf -uO $SMALL_NMSRT_BAM2.bam tmp-prefix2 | samtools bam2fq -s se2.fq.gz - | bwa mem -p /mnt/data/pipeline_genome_data/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta - > tmp.fq2
samtools bamshuf -uO $SMALL_NMSRT_BAM1.bam tmp-prefix | samtools bam2fq -s se.fq.gz - > tmp.fq1
samtools bamshuf -uO $SMALL_NMSRT_BAM2.bam tmp-prefix2 | samtools bam2fq -s se2.fq.gz - > tmp.fq2
rm -f tmp-prefix* se*.fq.gz
#samtools bam2fq $SMALL_NMSRT_BAM1.bam > tmp.fq1
#samtools bam2fq $SMALL_NMSRT_BAM2.bam > tmp.fq2
cat tmp.fq1 | grep '^@.*/1$' -A 3 --no-group-separator > $SMALL_FQ1_R1
cat tmp.fq1 | grep '^@.*/2$' -A 3 --no-group-separator > $SMALL_FQ1_R2
cat tmp.fq2 | grep '^@.*/1$' -A 3 --no-group-separator > $SMALL_FQ2_R1
cat tmp.fq2 | grep '^@.*/2$' -A 3 --no-group-separator > $SMALL_FQ2_R2
#sed -i 's/\/1/ 1:N:0:GATCAG/;N' $SMALL_FQ1_R1
#sed -i 's/\/2/ 2:N:0:GATCAG/;N' $SMALL_FQ1_R2
#sed -i 's/\/1/ 1:N:0:GATCAG/;N' $SMALL_FQ2_R1
#sed -i 's/\/2/ 2:N:0:GATCAG/;N' $SMALL_FQ2_R2
gzip -f $SMALL_FQ1_R1 $SMALL_FQ1_R2 $SMALL_FQ2_R1 $SMALL_FQ2_R2
rm -f $SMALL_BAM1 $SMALL_BAM2 $SMALL_NMSRT_BAM1.bam $SMALL_NMSRT_BAM2.bam tmp.fq1 tmp.fq2

BAM1=/srv/scratch/shared/surya/leepc12/run/NatProtocol_TestCase/ENCSR936XTK/PE/align/ctl1/ctl1-R1.PE2SE.bam
BAM2=/srv/scratch/shared/surya/leepc12/run/NatProtocol_TestCase/ENCSR936XTK/PE/align/ctl2/ctl2-R1.PE2SE.bam
SMALL_BAM1=$SAMPLE.$CHR.ctl1.bam
SMALL_BAM2=$SAMPLE.$CHR.ctl2.bam
SMALL_NMSRT_BAM1=$SAMPLE.$CHR.ctl1.nmsrt.bam
SMALL_NMSRT_BAM2=$SAMPLE.$CHR.ctl2.nmsrt.bam
SMALL_FQ1_R1=$SAMPLE.$CHR.ctl1-R1.fastq
SMALL_FQ1_R2=$SAMPLE.$CHR.ctl1-R2.fastq
SMALL_FQ2_R1=$SAMPLE.$CHR.ctl2-R1.fastq
SMALL_FQ2_R2=$SAMPLE.$CHR.ctl2-R2.fastq
mkdir -p $WORKDIR && cd $WORKDIR
samtools view -b $BAM1 "$CHR" > $SMALL_BAM1
samtools view -b $BAM2 "$CHR" > $SMALL_BAM2
samtools sort -n $SMALL_BAM1 $SMALL_NMSRT_BAM1
samtools sort -n $SMALL_BAM2 $SMALL_NMSRT_BAM2
samtools bamshuf -uO $SMALL_NMSRT_BAM1.bam tmp-prefix | samtools bam2fq -s se.fq.gz - > tmp.fq1
samtools bamshuf -uO $SMALL_NMSRT_BAM2.bam tmp-prefix2 | samtools bam2fq -s se2.fq.gz - > tmp.fq2
rm -f tmp-prefix* se*.fq.gz
#samtools bam2fq $SMALL_NMSRT_BAM1.bam > tmp.fq1
#samtools bam2fq $SMALL_NMSRT_BAM2.bam > tmp.fq2
cat tmp.fq1 | grep '^@.*/1$' -A 3 --no-group-separator > $SMALL_FQ1_R1
cat tmp.fq1 | grep '^@.*/2$' -A 3 --no-group-separator > $SMALL_FQ1_R2
cat tmp.fq2 | grep '^@.*/1$' -A 3 --no-group-separator > $SMALL_FQ2_R1
cat tmp.fq2 | grep '^@.*/2$' -A 3 --no-group-separator > $SMALL_FQ2_R2
#sed -i 's/\/1/ 1:N:0:TTAGGC/;N' $SMALL_FQ1_R1
#sed -i 's/\/2/ 2:N:0:TTAGGC/;N' $SMALL_FQ1_R2
#sed -i 's/\/1/ 1:N:0:TTAGGC/;N' $SMALL_FQ2_R1
#sed -i 's/\/2/ 2:N:0:TTAGGC/;N' $SMALL_FQ2_R2
gzip -f $SMALL_FQ1_R1 $SMALL_FQ1_R2 $SMALL_FQ2_R1 $SMALL_FQ2_R2
rm -f $SMALL_BAM1 $SMALL_BAM2 $SMALL_NMSRT_BAM1.bam $SMALL_NMSRT_BAM2.bam tmp.fq1 tmp.fq2
