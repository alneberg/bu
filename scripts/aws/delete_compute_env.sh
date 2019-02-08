VERSION=16
REGION="eu-west-1"
COMPUTE_ENV=Sarek_pricing_benchmark_${VERSION}
KEYNAME=Johannes_AWS
ACCOUNTID=$AWS_ACCOUNT_ID

# Disable and delete job queue
aws batch update-job-queue --region $REGION --job-queue ${COMPUTE_ENV} --state DISABLED
sleep 5
aws batch delete-job-queue --region $REGION --job-queue ${COMPUTE_ENV}
sleep 5

# Disable and delete compute env
aws batch update-compute-environment --region $REGION --compute-environment ${COMPUTE_ENV} --state DISABLED
sleep 5
aws batch delete-compute-environment --region $REGION --compute-environment ${COMPUTE_ENV}
