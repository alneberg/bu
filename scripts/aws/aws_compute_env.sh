VERSION=15
ACCOUNTID=$AWS_ACCOUNT_ID
SERVICEROLE=AWSBatchServiceRole
IAMFLEETROLE=AmazonEC2SpotFleetRole
SUBNETS=subnet-060e7d16f46ef3bd2,subnet-04cbdc281a40159ce,subnet-09c7a62f5f47124f0
SECGROUPS=sg-0e0d5e87597682f99
SPOTPER=50 # percentage of on demand
IMAGEID=ami-018adc348cc2bad99
INSTANCEROLE=ecsInstanceRole
KEYNAME=Johannes_AWS
MAXCPU=256 # max vCPUs in compute environment
INSTANCETYPES="i3"
COMPUTE_ENV=Sarek_pricing_benchmark_${VERSION}
TAGS={sarek_benchmark_compute=${VERSION},sarek_benchmark=${VERSION}}
REGION="eu-west-1"

# Create regular compute env and queue
aws batch create-compute-environment --region $REGION --compute-environment-name $COMPUTE_ENV --type MANAGED --state ENABLED --service-role ${SERVICEROLE} --compute-resources type=SPOT,minvCpus=0,maxvCpus=$MAXCPU,desiredvCpus=0,instanceTypes=$INSTANCETYPES,imageId=$IMAGEID,subnets=$SUBNETS,securityGroupIds=$SECGROUPS,ec2KeyPair=$KEYNAME,instanceRole=$INSTANCEROLE,bidPercentage=$SPOTPER,spotIamFleetRole=$IAMFLEETROLE,tags=$TAGS

sleep 5

aws batch create-job-queue --region $REGION --job-queue-name $COMPUTE_ENV --state ENABLED --priority 1 --compute-environment-order "order=1,computeEnvironment=$COMPUTE_ENV"

# Create tiny compute env and queue
# INSTANCETYPES="r4.large"
# COMPUTE_ENV=${COMPUTE_ENV}_tiny
 
# aws batch create-compute-environment --region $REGION --compute-environment-name $COMPUTE_ENV --type MANAGED --state ENABLED --service-role ${SERVICEROLE} --compute-resources type=SPOT,minvCpus=0,maxvCpus=$MAXCPU,desiredvCpus=0,instanceTypes=$INSTANCETYPES,imageId=$IMAGEID,subnets=$SUBNETS,securityGroupIds=$SECGROUPS,ec2KeyPair=$KEYNAME,instanceRole=$INSTANCEROLE,bidPercentage=$SPOTPER,spotIamFleetRole=$IAMFLEETROLE,tags=$TAGS

# sleep 5

# aws batch create-job-queue --region $REGION --job-queue-name $COMPUTE_ENV --state ENABLED --priority 1 --compute-environment-order "order=1,computeEnvironment=$COMPUTE_ENV"
