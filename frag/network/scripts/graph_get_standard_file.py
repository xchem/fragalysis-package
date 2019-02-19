#!/usr/bin/env python3
# coding=utf-8

"""A utility to get a standard data file from AWS S3.

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

parser = argparse.ArgumentParser('Graph Standard File Getter')
parser.add_argument('path', metavar='PATH', type=str,
                    help='The path, relative to the "standard" directory'
                         ' in your S3 bucket. e.g. "activity/senp7/standard-2"')
parser.add_argument('destination', metavar='DIR', type=str,
                    help='The local destination directory for the standard file,'
                         ' which must not exist and will be created')
parser.add_argument("--force", dest="force", action="store_true")

args = parser.parse_args()

# Destination?
if not args.force and os.path.isdir(args.destination):
    logger.error('Destination exists')
    sys.exit(1)

src = s3_standard_root + '/' + args.path + '/' + s3_standard_file
logger.info('Getting standard file from "%s"...', args.path)

# Get the contents of the selected directory...
s3_client = boto3.client('s3')
resp = s3_client.list_objects_v2(Bucket=s3_archive_bucket,
                                 Prefix=src)

# Download (if present)...
# An error if it isn't.
if resp and 'KeyCount' in resp and resp['KeyCount'] == 1:

    os.mkdir(args.destination)
    s3_client.download_file(s3_archive_bucket,
                            src,
                            os.path.join(args.destination, s3_standard_file))
    logger.info('Done')

else:

    logger.error('Standard file does not exist')
    sys.exit(1)
