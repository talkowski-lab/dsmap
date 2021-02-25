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
import tempfile
import argparse


# Set priority of dockers to build (in order)
# (Note: do not alter this without first confirming dependencies of base images)
ordered_dockers = ['athena', 'athena-cloud', 'dsmap', 'dsmap-cromwell']


class DockerError(Exception):
    """
    Dummy exception for errors in docker system calls
    """
    pass


class GithubError(Exception):
    """
    Dummy exception for errors in github system calls
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
    push_args_group.add_argument('--no-push', action='store_true', 
                                 help='Do not push images to GCR.')
    push_args_group.add_argument('--gcr-project', type=str, 
                                 help='GCR billing project to push the images to.')
    push_args_group.add_argument('--update-latest', action='store_true',
                                 help='also update \"latest\" tag in remote docker repo(s)')

    # Local args
    local_args_group = parser.add_argument_group('Local options')
    local_args_group.add_argument('--tmpdir', type=str, 
                                  help='Specify temporary directory to use as ' +
                                  'build context')

    return parser.parse_args()


def make_tmpdir(tmpdir_path):
    """
    Create temporary directory used as build context
    """

    if tmpdir_path is None:
        tmpdir_path = os.popen('mktemp -d').read().rstrip()
    else:
        os.mkdir(tmpdir_path)

    os.chdir(tmpdir_path)

    return tmpdir_path


def clone_git_repo(repo_url, hash):
    """
    Clone a repo from GitHub and checkout specified hash
    """

    repo_name = os.path.basename(repo_url).split('.')[0]

    # Clone
    ret = os.system('git clone ' + repo_url)
    if 0 != ret:
        raise GithubError('Failed to clone ' + repo_url)

    # Checkout
    checkout_cmd = 'cd ' + repo_name + ' && '
    checkout_cmd += 'git checkout ' + hash + ' && '
    checkout_cmd += 'cd -'
    ret = os.system(checkout_cmd)
    if 0 != ret:
        raise GithubError('Failed to checkout ' + hash)


def format_build_args(args, docker, build_info, dockers_to_build):
    """
    Extract & format image-specific docker build arguments
    """

    build_args = ' '

    if docker == 'athena-cloud':
        # Use most recent build of athena base image if athena was also rebuilt
        if 'athena' in dockers_to_build:
            build_args += '--build-arg ATHENA_BASE_IMAGE={} '.format(build_info['athena']['remote'])

    elif docker == 'dsmap':
        # Use most recent build of athena-cloud base image if athena-cloud was also rebuilt
        if 'athena-cloud' in dockers_to_build:
            build_args += '--build-arg ATHENA_CLOUD_BASE_IMAGE={} '.format(build_info['athena-cloud']['remote'])

    elif docker == 'dsmap-cromwell':
        # Use most recent build of dsmap base image if dsmap was also rebuilt
        if 'dsmap' in dockers_to_build:
            build_args += '--build-arg DSMAP_BASE_IMAGE={} '.format(build_info['dsmap']['remote'])

    build_args += '--tag {} '.format(build_info[docker]['remote'])

    return build_args


def build_docker(build_info, build_args):
    """
    Build a single Docker image
    """

    remote = build_info['remote']
    dockerfile_dir = build_info['dockerfile_dir']

    # Construct docker build call
    build_cmd = 'docker build --progress plain ' + build_args 
    build_cmd += '-f {}/Dockerfile '.format(dockerfile_dir)
    build_cmd += tmpdir_path + '\n'

    # Build docker
    print('\nNow building ' + remote + ' as follows:\n')
    print(build_cmd)
    ret = os.system(build_cmd)
    if 0 != ret:
        raise DockerError('Failed to build ' + remote)
    else:
        print('\nSuccessfully built {}\n'.format(remote))


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
        print('\nSuccessfully pushed {}\n'.format(remote))


def retag_docker(old_remote, new_remote):
    """
    Tags a built docker to a new remote
    """

    # Construct docker tag call
    tag_cmd = 'docker tag {} {}'.format(old_remote, new_remote)

    # Retag image
    ret = os.system(tag_cmd)
    if 0 != ret:
        raise DockerError("Failed to retag {} as {}".format(old_remote, new_remote))
    else:
        print('\nSuccessfully retagged {} as {}\n'.format(old_remote, new_remote))


def cleanup():
    """
    Clean up build context and switch back to execution directory
    """

    os.chdir(exec_dir)
    os.system('rm -rf ' + tmpdir_path)


def main():
    # Read command-line arguments
    args = parse_commandline_args()

    # Get path to local DSMap repo
    global exec_dir
    exec_dir = os.getcwd()
    dsmap_dir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-2])

    # Determine which images to build
    if 'All' in args.images:
        dockers_to_build = ordered_dockers
    else:
        dockers_to_build = [a for a in ordered_dockers if a in args.images]

    # Resolve paths to Dockerfiles
    print('\nBuilding the following Docker images:\n')
    for docker in dockers_to_build:
        print('  - {}:{}\n'.format(docker, args.tag))
    build_info = {d : {'dockerfile_dir' : '{}/dockerfiles/{}'.format(dsmap_dir, d),
                       'remote' : 'us.gcr.io/{}/{}:{}'.format(args.gcr_project, d, args.tag)} \
                  for d in dockers_to_build}

    # Create temporary directory to use as build context
    try:
        global tmpdir_path
        tmpdir_path = make_tmpdir(args.tmpdir)
        print('\nUsing {} as temporary build context\n'.format(tmpdir_path))

        # Clone necessary repos
        print('\nCloning required GitHub repos into local build context\n')
        if 'athena' in dockers_to_build:
            clone_git_repo('git@github.com:talkowski-lab/athena.git', args.athena_hash)
        if 'dsmap' in dockers_to_build:
            clone_git_repo('git@github.com:talkowski-lab/dsmap.git', args.dsmap_hash)

        # Build each Docker
        for docker in dockers_to_build:
            build_args = format_build_args(args, docker, build_info, dockers_to_build)
            build_docker(build_info[docker], build_args)

        # Push each built Docker if GCR project is provided
        if args.gcr_project is not None and not args.no_push:
            print('\nPushing the following Docker images:\n')
            for docker in dockers_to_build:
                print('  - {}\n'.format(build_info[docker]['remote']))

            for docker in dockers_to_build:
                push_docker(build_info[docker]['remote'])

            # Also update :latest tag, if optioned
            if args.update_latest and args.tag != 'latest':
                print('\nAlso updating latest tag for the following Docker images:\n')
                for docker in dockers_to_build:
                    print('  - {} => us.gcr.io/{}/{}:latest\n'.format(build_info[docker]['remote'], args.gcr_project, docker))
                for docker in dockers_to_build:
                    newtag = 'us.gcr.io/{}/{}:latest'.format(args.gcr_project, docker)
                    retag_docker(build_info[docker]['remote'], newtag)
                    push_docker(newtag)

        # Clean up, if successful
        cleanup()
        print('\nFinished building all docker images\n')

    except:
        cleanup()


if __name__ == '__main__':
    main()
