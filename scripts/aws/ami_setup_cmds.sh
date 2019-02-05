sudo yum -y update

# These didn't seem to be necessary...
sudo systemctl stop ecs
sudo rm /var/lib/ecs/data/ecs_agent_data.json


# Install conda
sudo yum install -y bzip2 wget
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $HOME/miniconda
$HOME/miniconda/bin/conda install -c conda-forge -y awscli
rm Miniconda3-latest-Linux-x86_64.sh

# Setup script to be ran on boot
sudo yum install -y mdadm # Already installed
wget https://raw.githubusercontent.com/alneberg/bu/master/scripts/aws/mount_physical_drives_aws.py
sudo cp mount_physical_drives_aws.py /usr/local/bin/

sudo sh -c "echo -e '[Unit]\nDescription=Convert free block devices into raid and add to Docker\nBefore=ecs.service\n[Service]\nExecStart=/usr/local/bin/mount_physical_drives_aws.py\n[Install]\nWantedBy=multi-user.target' > /etc/systemd/system/mount_drives.service"

sudo systemctl enable mount_drives
