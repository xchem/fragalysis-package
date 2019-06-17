#!/usr/bin/env python3
# coding=utf-8

"""A utility to write a build provenance files to AWS S3.
All files in the source directory that end '.prov' will be written.

To use this utility you will need to ensure that your AWS credentials
are available via expected environment variables. Please refer to
the AWS boto3 documentation, which describes this: -

    https://boto3.amazonaws.com/v1/documentation/api/latest/guide/quickstart.html

Alan Christie
March 2019
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
s3_standard_root = 'standard'
s3_archive_bucket = os.environ.get(s3_bucket_env)
if not s3_archive_bucket:
    logger.error('You must define %s', s3_bucket_env)
    sys.exit(1)

parser = argparse.ArgumentParser('Graph Standard Provenance File Putter')
parser.add_argument('source', metavar='DIR', type=str,
                    help='The local directory (where the .prov files exist)')
parser.add_argument('path', metavar='PATH', type=str,
                    help='The path, relative to the "standard" directory'
                         ' in your S3 bucket of the standard.'
                         'e.g. "activity/senp7/standard-2"')

args = parser.parse_args()

# Check source directory
if not os.path.isdir(args.source):
    logger.error('The source directory does not exist (%s)', args.source)
    sys.exit(1)

s3_client = boto3.client('s3')

# Upload the provenance files in the source directory...
potential_files = os.listdir(args.source)
for potential_file in potential_files:
    src = os.path.join(args.source, potential_file)
    if os.path.isfile(src) and potential_file.endswith('.prov'):
        dst = s3_standard_root + '/' + args.path + '/' + potential_file
        logger.info('Putting %s -> %s...', potential_file, args.path)
        s3_client.upload_file(src, s3_archive_bucket, dst)
