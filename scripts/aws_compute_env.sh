VERSION=9
ACCOUNTID=$AWS_ACCOUNT_ID
SERVICEROLE=AWSBatchServiceRole
IAMFLEETROLE=AmazonEC2SpotFleetRole
SUBNETS=subnet-060e7d16f46ef3bd2,subnet-04cbdc281a40159ce,subnet-09c7a62f5f47124f0
SECGROUPS=sg-0e0d5e87597682f99
SPOTPER=50 # percentage of on demand
IMAGEID=ami-0978cc20a5c1b017e
INSTANCEROLE=ecsInstanceRole
KEYNAME=Johannes_AWS
MAXCPU=256 # max vCPUs in compute environment
INSTANCETYPES="m4.4xlarge,m4.xlarge,r4.4xlarge,r4.2xlarge,r4.xlarge,m5.4xlarge"
COMPUTE_ENV=Sarek_pricing_benchmark_${VERSION}
TAGS={Sarek-benchmark-compute=${VERSION},Sarek-benchmark=${VERSION}}

aws batch create-compute-environment --compute-environment-name $COMPUTE_ENV --type MANAGED --state ENABLED
--service-role ${SERVICEROLE} --compute-resources type=SPOT,minvCpus=0,maxvCpus=$MAXCPU,desiredvCpus=0,instanceTypes=$INSTANCETYPES,imageId=$IMAGEID,subnets=$SUBNETS,securityGroupIds=$SECGROUPS,ec2KeyPair=$KEYNAME,instanceRole=$INSTANCEROLE,bidPercentage=$SPOTPER,spotIamFleetRole=$IAMFLEETROLE,tags=$TAGS

sleep 5

aws batch create-job-queue --job-queue-name $COMPUTE_ENV --state ENABLED --priority 1 --compute-environment-order "order=1,computeEnvironment=$COMPUTE_ENV"
