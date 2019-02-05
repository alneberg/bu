#!/usr/bin/env python
"""A script to transfer git repositories to secure uppmax cluster irma

Uses the ssh config file ~/.ssh/config and an additional ~/.remote_hosts.yaml config file.
"""
import argparse
import os
from os.path import join as ospj
import sys
from fabric import Connection
from git import Repo
import coloredlogs
import logging
import yaml
from yaml import Loader

# Set up a logger with colored output
logger = logging.getLogger('transfer_repo_to_irma')
logger.propagate = False  # Otherwise the messages appeared twice
coloredlogs.install(level='INFO', logger=logger,
                            fmt='%(asctime)s %(levelname)s %(message)s')

def get_dirname(outdir):
    if outdir == '.':
        outdir = os.getcwd()
    basename = os.path.basename(os.path.normpath(outdir))
    return basename

def main(args):
    with open(os.path.expanduser(args.config_file), 'r') as fh:
        config = yaml.load(fh)

    # Create an archive
    dirname = get_dirname(args.repo)
    
    logger.info("Using repo {}".format(dirname))

    repo = Repo(args.repo)
    assert not repo.is_dirty()

    archive_name = dirname
    git_tag = next((tag for tag in repo.tags if tag.commit == repo.head.commit), None)
    if git_tag:
        archive_name += '_' + git_tag    

    else:
        archive_name += '_' + repo.head.object.hexsha

    if args.extra_tag:
        archive_name += '_' + args.extra_tag

    logger.info("Creating repo archive {}".format(archive_name))
    archive = "{0}.tar.gz".format(archive_name)
    archive_path = ospj(args.repo, archive)
    with open(archive_path, 'wb') as fh:
        repo.archive(fh, format='tar.gz')

    logger.info("Archive created.")

    # Transfer archive to remote
    remote_dir = config['hosts']['irma']['archive_dir']

    Connection('irma').put(archive_path, remote=remote_dir)
    logger.info("Archive successfully transferred to irma")

    # Extract remote archive
    c = Connection('irma')     
    remote_archive_path = ospj(remote_dir, archive)
    remote_extracted_path = remote_archive_path.replace('.tar.gz', '')
    c.run('mkdir {}'.format(remote_extracted_path))
    c.run('tar -xvzf {} -C {}'.format(remote_archive_path, remote_extracted_path))

    # Create a link from dev or latest to the new archive
    if args.mode == 'dev': 
        link_name = "{}_dev".format(dirname)
    else:
        link_name = "{}_latest".format(dirname)
    c.run('cd {}; ln -s {} {}'.format(remote_dir, remote_extracted_path, link_name))
    logger.info("{} successfully linked as the new {}".format(dirname, link_name))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', default='.', help=("Repository that will be archived and "
            "transfered to remote. Default is current directory."))
    parser.add_argument('--extra_tag', default="", help=("Extra tag to be appended to the "
            "archive name in case the latest git tag is not sufficient"))
    parser.add_argument('--config_file', default='~/.remote_hosts.yaml', help=('Config file '
                'for file paths on remote server, default ~/.remote_hosts.yaml'))
    parser.add_argument('--mode', default='dev', choices=['dev', 'prod'], help=('Choose mode of '
        'publication: dev if to {repo}_dev and prod if to {repo}_latest. Default is dev.'))
    args = parser.parse_args()
    main(args)
