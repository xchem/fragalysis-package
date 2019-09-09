#!/usr/bin/env python3
# coding=utf-8

"""A utility to write a standard file to AWS S3.

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
s3_bucket_env = 'AWS_S3_BUCKET'
s3_standard_root = 'standard'
s3_standard_file = 'standardised-compounds.tab.gz'
s3_archive_bucket = os.environ.get(s3_bucket_env)
if not s3_archive_bucket:
    logger.error('You must define %s', s3_bucket_env)
    sys.exit(1)

parser = argparse.ArgumentParser('Graph Standard File Putter')
parser.add_argument('source', metavar='DIR', type=str,
                    help='The local directory (where the standard file exists)')
parser.add_argument('path', metavar='PATH', type=str,
                    help='The path, relative to the "standard" directory'
                         ' in your S3 bucket. e.g. "activity/senp7/standard-1"')

args = parser.parse_args()

# Standard file
filename = os.path.join(args.source, s3_standard_file)
if not os.path.isfile(filename):
    logger.error('Your standard file (%s) does not exist', filename)
    sys.exit(1)

s3_client = boto3.client('s3')

# We must not write to an existing path (key).
# It is up to the user to make sure the destination does not exist,
# it's too easy to over-write files in S3.

logger.info('Putting %s -> %s...', s3_standard_file, args.path)

dst = s3_standard_root + '/' + args.path + '/' + s3_standard_file
target = s3_client.list_objects_v2(Bucket=s3_archive_bucket,
                                   Prefix=dst)
if 'KeyCount' in target and target['KeyCount']:
    logger.error('The standard file already exists in S3.'
                 ' You cannot "put" to existing locations.')
    sys.exit(1)

# Upload the standard file...
s3_client.upload_file(filename, s3_archive_bucket, dst)
