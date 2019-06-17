#!/usr/bin/env python3
# coding=utf-8

"""A utility to write graph files to AWS S3.

Graph files consist of everything (except provenance files) in the names
directory. The file are the results os running the 'process' scripts.
The files are deposited into the build directory

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
s3_bucket_env = 'FRAGALYSIS_S3_BUCKET'
s3_build_root = 'build'
s3_archive_bucket = os.environ.get(s3_bucket_env)
if not s3_archive_bucket:
    logger.error('You must define %s', s3_bucket_env)
    sys.exit(1)

parser = argparse.ArgumentParser('Graph Build Graph File Putter')
parser.add_argument('source', metavar='DIR', type=str,
                    help='The local directory (where the graph data exists)')
parser.add_argument('path', metavar='PATH', type=str,
                    help='The path, relative to the "raw" directory'
                         ' in your S3 bucket. e.g. "activity/senp7/v1/build-1"')

args = parser.parse_args()

# Check source directory
if not os.path.isdir(args.source):
    logger.error('The source directory does not exist (%s)', args.source)
    sys.exit(1)

# The S3 build path...
dst = s3_build_root + '/' + args.path + '/'
s3_client = boto3.client('s3')

# We write everything in the source directory that's a file...
num_files_saved = 0
items = os.listdir(args.source)
for item in items:
    # Skip any .prov files.
    if item.endswith('.prov'):
        continue
    src = os.path.join(args.source, item)
    if os.path.isfile(src):
        dst = s3_build_root + '/' + args.path + '/' + item
        logger.info('Putting %s -> %s...', item, args.path)
        s3_client.upload_file(src, s3_archive_bucket, dst)
        num_files_saved += 1

# It's an error not to have saved any files
if num_files_saved == 0:
    logger.error('No files were saved, is the directory empty?')
    sys.exit(2)
