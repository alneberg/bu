#!/bin/env python
"""Convert unmounted and unpartioned disks to raid0

Following recipes from https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/raid-config.html

Requires mdadm to be installed
"""

import sys
import time
import os
import shutil

import json
import subprocess


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2019 Genome Institute of Singapore"
__license__ = "The MIT License (MIT)"


LABEL = "myraid"


def get_free_disks():
    """List unmounted and unpartitioned disks
    """
    res = subprocess.check_output(['lsblk', '-J'], stderr=subprocess.STDOUT)
    jd = json.loads(res.decode())
    devs = []
    for bd in jd['blockdevices']:
        if not bd['mountpoint'] and not 'children' in bd:
            devs.append("/dev/{}".format(bd["name"]))
    return devs



def create_raid0(disks):
    """main procedure
    """

    print("# Creating raid0 array")
    # use '--force' in case we only have 1 device
    cmd = ("mdadm --create --verbose /dev/md0 --force --level=0"
           " --name={} --raid-devices={} {}".format(LABEL, len(disks), ' '.join(disks))
          )
    res = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    print(res.decode())

    print("# Waiting for sync")
    # FIXME crutch. could try to repeatedly read /proc/mdstat until e.g resync disappears
    time.sleep(10)

    print("# Creating config file for reassembly on boot")
    cmd = "mdadm --detail --scan | sudo tee -a /etc/mdadm.conf"
    res = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    print(res.decode())

    print("# Creating ramdisk to preload modules")
    cmd = "dracut -H -f /boot/initramfs-$(uname -r).img $(uname -r)"
    res = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    print(res.decode())

    print("# Creating FS")
    cmd = "mkfs.ext4 -L {} /dev/md0".format(LABEL)
    res = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    print(res.decode())

    return "/dev/md0"


def dev_to_docker(dev):
    print("# Stopping Docker and clearing root")
    cmd = "systemctl stop docker"
    res = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    print(res.decode())
    shutil.rmtree('/var/lib/docker')
    os.mkdir('/var/lib/docker')

    print("# Creating FS")
    cmd = "mkfs.ext4 -L {} {}".format(LABEL, dev)
    res = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    print(res.decode())

    print("# Adding device to fstab and mounting")
    with open("/etc/fstab", 'a') as fh:
        fh.write("\n{}\t/var/lib/docker\text4\tdefaults,nofail\t0\t2\n".format(dev))
    cmd = "mount /var/lib/docker"
    res = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    print(res.decode())

    print("# Starting Docker")
    cmd = "systemctl start docker"
    res = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    print(res.decode())

    print("# Starting ECS")
    cmd = "systemctl start ecs"
    res = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    print(res.decode())

def main():
    """main procedure
    """

    disks = get_free_disks()
    if not disks:
        sys.stderr.write("No suitable devices found\n")
        sys.exit(1)
    print("# Free disks: {}".format(', '.join(disks)))

    raid = create_raid0(disks)

    dev_to_docker(raid)


if __name__ == "__main__":
    main()

