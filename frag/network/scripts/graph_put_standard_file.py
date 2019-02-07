#!/usr/bin/env python
# coding=utf-8

"""A utility to write standard file to AWS S3.
This utility saves data, normally after initial standardising.

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
s3_standard_file = 'standardised-compounds.tab.gz'
s3_archive_bucket = os.environ.get(s3_bucket_env)
if not s3_archive_bucket:
    logger.error('You must define %s', s3_bucket_env)
    sys.exit(1)

parser = argparse.ArgumentParser('Graph Standard Data Putter')
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

# The filename of the standard?
dst = s3_standard_root + '/' + args.path + '/' + s3_standard_file
logger.info('Putting "%s"...', dst)

# Upload the standard file...
s3_client = boto3.client('s3')
s3_client.upload_file(filename, s3_archive_bucket, dst)
