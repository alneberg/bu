#!/bin/bash

# Source: @marcelm: https://gist.github.com/marcelm/2aa3c9c1678f39e6793f8c000d72176c
#
# A workaround for an issue with Nextflow (which may actually be a bash bug),
# see <https://github.com/SciLifeLab/Sarek/issues/420>
#
# The problem is that Nextflow does not notice that a job has finished and
# hangs indefinitely.
#
# This script looks for zombie processes that are children of a script named
# .command.stub, and kills that script. This seems to let the pipeline continue
# in a normal way. No guarantees that there are no side-effects.
#
# To run this for all running jobs of a certain user on a SLURM cluster,
# something like this can be used (kill-zombie.sh needs to be in the $PATH, or
# use the absolute path name):
#
# squeue -h -t R -o %B -u $USER | sort -u | xargs -n1 -P20 sh -c 'ssh $1 kill-zombie.sh' --

set -euo pipefail
pgrep java > /dev/null || exit

# PIDs of all zombies
zombie_pids=($(ps -eo stat=,pid=|grep '^Z'|awk '{print $2}'))

# Wait a while to give "legitimate" zombies time to disappear
sleep 1m

# Iterate over all zombie PIDs
for zombie in ${zombie_pids[*]}; do
    # check whether the zombie is still there
    if ps -q $zombie -o stat=|grep -q Z; then
        # Find out the parent's PID
        parent=$(ps -eo ppid= -q $zombie)
        # Kill the parent only if it is named '.command.stub'
        if ps -p $parent -o args=|grep -q -F .command.stub; then
            # this is a good candidate
            echo "Killing parent of zombie on $HOSTNAME"
            kill -9 $parent
        fi
    fi
done
