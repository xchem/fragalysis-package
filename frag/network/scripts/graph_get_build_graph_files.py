#!/usr/bin/env python3
# coding=utf-8

"""A utility to get all the graph files from a build on AWS S3.

This module assumes that the source of the data
resides in an S3 bucket with a directory structure that complies with the
**Fragalysis S3 Data Archive** definition.

To use this utility you will need to ensure that your AWS credentials
are available via expected environment variables. Please refer to
the AWS boto3 documentation, which describes this: -

    https://boto3.amazonaws.com/v1/documentation/api/latest/guide/quickstart.html

Alan Christie
February 2019
"""

import argparse
import logging
import os
import re
import stat
import sys

import boto3

# Configure basic logging
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')

out_hdlr = logging.StreamHandler()
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(out_hdlr)

# Expected environment variables (that define the bucket)
s3_bucket_env = 'AWS_S3_BUCKET'
s3_build_root = 'build'
s3_archive_bucket = os.environ.get(s3_bucket_env)
if not s3_archive_bucket:
    logger.error('You must define %s', s3_bucket_env)
    sys.exit(1)

parser = argparse.ArgumentParser('Graph Graph Build File Getter')
parser.add_argument('path', metavar='PATH', type=str,
                    help='The path, relative to the "build" directory'
                         ' in your S3 bucket. e.g. "activity/senp7/v1/build-1"')
parser.add_argument('destination', metavar='DIR', type=str,
                    help='The local destination directory for the data,'
                         ' which must not exist and will be created')
parser.add_argument("--force", action="store_true",
                    help='Force retrieval of the files, even if the'
                         ' destination directory exists')
parser.add_argument("--for-combination", action="store_true",
                    help='Only get the files required for combination'
                         ' (i.e. the de-duplicated node and edge CSV files')

args = parser.parse_args()

# Destination?
if not args.force and os.path.isdir(args.destination):
    logger.error('Destination exists')
    sys.exit(1)

src = s3_build_root + '/' + args.path + '/'
logger.info('Getting graph files from "%s"...', args.path)

# Get the contents of the selected directory...
s3_client = boto3.client('s3')
resp = s3_client.list_objects_v2(Bucket=s3_archive_bucket,
                                 Prefix=src)

# A reg-ex list of files to exclude
exclude_files = [r'nodes.*',
                 r'.*[.]prov',
                 r'.*[.]txt',
                 r'excluded[.].*',
                 r'.*[.]txt[.]gz',
                 r'done']
combination_files = ['edges.csv.gz',
                     'nodes.csv.gz']

# And download everything that looks like a file...
if resp and 'KeyCount' in resp:

    if resp['KeyCount']:

        if not os.path.isdir(args.destination):
            os.mkdir(args.destination)

        num_files = 0
        for content in resp['Contents']:
            if not content['Key'].endswith('/'):
                filename = content['Key'].split('/')[-1]
                # Is the file excluded?
                get_file = True
                if not args.for_combination:
                    for regex in exclude_files:
                        if re.compile(regex).match(filename):
                            get_file = False
                            break
                else:
                    if filename not in combination_files:
                        get_file = False
                # Get the file?
                if get_file:
                    logger.info('%s > %s', filename, args.destination)
                    s3_client.download_file(s3_archive_bucket,
                                            content['Key'],
                                            os.path.join(args.destination, filename))
                    num_files += 1
        logger.info('Done (%d)', num_files)

        # Change permission of .py and .sh files...
        downloaded_files = os.listdir(args.destination)
        for downloaded_file in downloaded_files:
            if downloaded_file.endswith('.sh') or downloaded_file.endswith('.py'):
                permissions = stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR
                permissions |= stat.S_IRGRP | stat.S_IXGRP
                permissions |= stat.S_IROTH | stat.S_IXOTH
                os.chmod(os.path.join(args.destination, downloaded_file), permissions)

    else:

        # KeyCount was zero
        logger.warning('No files found')
        sys.exit(1)

else:

    logger.error('Unexpected response (expected KeyCount)')
    sys.exit(1)
