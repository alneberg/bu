# Based on https://gist.github.com/andreas-wilm/3460a788d6548370a136e63b5b91281e
# by @andreas-wilm

# This script automates the mounting of AWS attached storage into raid0.
# https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/raid-config.html

# Not including xvda as this is the root partition
nvmes=$(sudo lsblk | awk '/^nvme/ {printf "/dev/%s ", $1} ')
sudo mdadm --create --verbose /dev/md0 --level=0 --name=my_raid --raid-devices=$(echo $nvmes | wc -w) $nvmes
sleep 10; 
# Check if the file system exists (the instance might have been up for some time already)
fs_check=`sudo file -s /dev/md0`
if [[ $fs_check = "/dev/md0: data" ]]; then
    sudo mkfs.ext4 -L my_raid /dev/md0
    sudo mdadm --detail --scan | sudo tee -a /etc/mdadm.conf
    sudo dracut -H -f /boot/initramfs-$(uname -r).img $(uname -r)
else
    # In case this job was started at the same time as the first one
    sleep 10;
fi
sudo mkdir /data
sudo mount LABEL=my_raid /data
sudo chown ec2-user:ec2-user /data/
