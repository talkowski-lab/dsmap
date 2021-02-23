#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2021-Present Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Helper tool to rebuild Athena and DSMap Docker images
"""

# Note: loosely based on GATK-SV build_docker.py


import sys, argparse, os, os.path
import argparse


# Set priority of dockers to build (in order)
# (Note: do not alter this without first confirming dependencies of base images)
ordered_dockers = ['athena', 'athena-cloud', 'dsmap']


class DockerError(Exception):
    """
    Dummy exception for errors in docker system calls
    """
    pass


def parse_commandline_args():
    """
    Wrapper for handling command-line arguments
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Required args
    required_args_group = parser.add_argument_group('Required')

    # Build args
    build_args_group = parser.add_argument_group('Build options', 'Options when building images')
    required_args_group.add_argument('-t', '--tag', type=str, default="latest",
                                     help="Tag applied to all images.")
    build_args_group.add_argument('-i', '--images', nargs='+', type=str, default='All',
                                  choices=['All'] + ordered_dockers, 
                                  help='Specify which images to build [default: All]')

    # GitHub args
    github_args_group = parser.add_argument_group('GitHub options', 'Options for GitHub repos')
    required_args_group.add_argument('-a', '--athena-hash', type=str, default="master",
                                     help="GitHub commit hash for version of Athena to build.")
    required_args_group.add_argument('-d', '--dsmap-hash', type=str, default="main",
                                     help="GitHub commit hash for version of DSMap to include.")

    # Docker push args
    push_args_group = parser.add_argument_group('Docker push options',
                                                'Options to push images to remote container registry')
    push_args_group.add_argument('--gcr-project', type=str, 
                                 help='GCR billing project to push the images to.')
    push_args_group.add_argument('--update-latest', action='store_true',
                                 help='also update \"latest\" tag in remote docker repo(s)')

    return parser.parse_args()


def format_build_args(args, docker, build_info):
    """
    Extract & format image-specific docker build arguments
    """

    build_args = ' '

    if docker == 'athena':
        build_args += '--build-arg ATHENA_COMMIT={} '.format(args.athena_hash)
    elif docker == 'dsmap':
        build_args += '--build-arg ATHENA_BASE_IMAGE=us.gcr.io/broad-dsmap/athena:{} '.format(args.athena_hash)
        build_args += 'DSMAP_COMMIT={} '.format(args.dsmap_hash)

    build_args += '--tag {} '.format(build_info['remote'])

    return build_args


def build_docker(build_info, build_args):
    """
    Build a single Docker image
    """

    remote = build_info['remote']
    dockerfile_dir = build_info['dockerfile_dir']

    # Construct docker build call
    build_cmd = 'cd {} && '.format(dockerfile_dir)
    build_cmd += 'docker build --progress plain ' + build_args + ' . && cd - \n'

    # Build docker
    print('Now building ' + remote + ' as follows:')
    print(build_cmd)
    ret = os.system(build_cmd)
    if 0 != ret:
        raise DockerError("Failed to build " + remote)
    else:
        print('Successfully built ' + remote)


def push_docker(remote):
    """
    Push a single tagged docker image
    """

    # Construct docker push call
    push_cmd = 'docker push ' + remote

    # Push image
    ret = os.system(push_cmd)
    if 0 != ret:
        raise DockerError("Failed to push " + remote)
    else:
        print('Successfully pushed ' + remote)


def retag_docker(old_remote, new_remote):
    """
    Tags a built docker to a new remote
    """

    # Construct docker tag call
    tag_cmd = 'docker tag {} {}'.format(old_remote, new_remote)

    # Retag image
    ret = os.system(push_cmd)
    if 0 != ret:
        raise DockerError("Failed to retag {} as {}".format(old_remote, new_remote))
    else:
        print('Successfully retagged {} as {}'.format(old_remote, new_remote))


def main():
    # Read command-line arguments
    args = parse_commandline_args()

    # Get path to local DSMap repo
    dsmap_dir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-2])

    # Resolve paths to dockerfiles
    dockers_to_build = [a for a in ordered_dockers if a in args.images]
    print('\nBuilding the following Docker images:\n')
    for docker in dockers_to_build:
        print('  - {}:{}\n'.format(docker, args.tag))
    build_info = {d : {'dockerfile_dir' : '{}/dockerfiles/{}'.format(dsmap_dir, d),
                       'remote' : 'us.gcr.io/{}/{}:{}'.format(args.gcr_project, d, args.tag)} \
                  for d in dockers_to_build}

    # Build each Docker
    for docker in dockers_to_build:
        build_args = format_build_args(args, docker, build_info[docker])
        build_docker(build_info[docker], build_args)

    # Push each built Docker if GCR project is provided
    if args.gcr_project is not None:
        print('\nPushing the following Docker images:\n')
        for docker in dockers_to_build:
            print('  - {}\n'.format(build_info[docker]['remote']))

        for docker in dockers_to_build:
            push_docker(build_info[docker]['remote'])

        # Also update :latest tag, if optioned
        if args.update_latest and args.tag != 'latest':
            print('\nAlso updating latest tag for the following Docker images:\n')
            newtag = 'us.gcr.io/{}/{}:latest'.format(args.gcr_project, docker)
            retag_docker(build_info['docker']['remote'], newtag)
            push_docker(newtag)


if __name__ == '__main__':
    main()
