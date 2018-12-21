WORKDIR="s3://sarek-work-benchmark/work-tiny_2018_12_13"
OUTDIR="s3://sarek-result-benchmark/results_tiny_20181213"
GENOME_BASE="s3://sarek-references/small/"
REPORT_DIR="../Reports_tiny"

COMMON_PARAMS="-profile awsbatch --awsqueue $AWS_QUEUE --awsqueue_tiny $AWS_QUEUE_TINY -work-dir $WORKDIR --outDir
$OUTDIR --verbose -resume"

STEP=main
echo $STEP
echo nextflow run main.nf $COMMON_PARAMS --localReportdir ${REPORT_DIR}_${STEP} --sample Sarek-data/testdata/tsv/tiny-s3.tsv --genome smallGRCh37 --genome_base $GENOME_BASE 

STEP=germlineVC
echo $STEP
echo nextflow run germlineVC.nf $COMMON_PARAMS --localReportdir ${REPORT_DIR}_${STEP} --sample $OUTDIR/Preprocessing/Recalibrated/recalibrated.tsv --genome smallGRCh37 --genome_base $GENOME_BASE --tools HaplotypeCaller,Manta,Strelka

STEP=annotate
echo $STEP
echo nextflow run annotate.nf $COMMON_PARAMS --localReportdir ${REPORT_DIR}_${STEP} --genome smallGRCh37 --genome_base $GENOME_BASE --annotateTools HaplotypeCaller,Strelka --tools snpEFF,VEP

STEP=multiqc
echo $STEP
echo nextflow run runMultiQC.nf $COMMON_PARAMS --localReportdir ${REPORT_DIR}_${STEP}

WORKDIR="s3://sarek-work-benchmark/work-dream-ssd_2018_12_12"
OUTDIR="s3://sarek-result-benchmark/results_dream_ssd_20181212"
GENOME_BASE="s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38"
AWS_QUEUE="Sarek_pricing_benchmark_10"
AWS_QUEUE_TINY="Sarek_pricing_benchmark_10_tiny"
REPORT_DIR="Reports"


