set -e
WORKDIR="s3://sarek-work-benchmark/work_190205"
OUTDIR="s3://sarek-result-benchmark/results_190205"
GENOME_BASE="s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38"
AWS_QUEUE="Sarek_pricing_benchmark_16"
AWS_QUEUE_TINY="Sarek_pricing_benchmark_16"
REPORT_DIR="../Reports_dream_m5d_190205"
MAIN_SAMPLE_TSV="Sarek-data/testdata/tsv/dream-test-normal-s3.tsv"
GENOME="GRCh38"

#Comment away this when running real data
#WORKDIR="s3://sarek-work-benchmark/work-tiny_20181221"
#OUTDIR="s3://sarek-result-benchmark/results_tiny_20181221"
#GENOME_BASE="s3://sarek-references/small/"
#REPORT_DIR="../Reports_tiny"
#MAIN_SAMPLE_TSV="Sarek-data/testdata/tsv/tiny-s3.tsv"
#GENOME="smallGRCh37"

COMMON_PARAMS="-profile awsbatch --awsqueue $AWS_QUEUE --awsqueue_tiny $AWS_QUEUE_TINY -work-dir $WORKDIR --outDir $OUTDIR --verbose -resume"


STEP=main
echo $STEP
CMD="nextflow run main.nf $COMMON_PARAMS --localReportDir ${REPORT_DIR}_${STEP} --sample $MAIN_SAMPLE_TSV  --genome $GENOME --genome_base $GENOME_BASE"
echo $CMD
$CMD

STEP=germlineVC
echo $STEP
nextflow run germlineVC.nf $COMMON_PARAMS --localReportDir ${REPORT_DIR}_${STEP} --sample $OUTDIR/Preprocessing/Recalibrated/recalibrated.tsv --genome $GENOME --genome_base $GENOME_BASE --tools HaplotypeCaller,Manta,Strelka

STEP=annotate
echo $STEP
nextflow run annotate.nf $COMMON_PARAMS --localReportDir ${REPORT_DIR}_${STEP} --genome $GENOME --genome_base $GENOME_BASE --annotateTools HaplotypeCaller,Strelka --tools snpEFF,VEP

STEP=multiqc
echo $STEP
nextflow run runMultiQC.nf $COMMON_PARAMS --localReportDir ${REPORT_DIR}_${STEP}

