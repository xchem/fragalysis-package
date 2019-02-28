#!/usr/bin/env python3
# coding=utf-8

"""A utility to check that a given build path is present.
Call with something like "build/activity/senp7/build-1". If it is present
the output of the utility is "present" if it is not present
the output is "absent".

This utility must be given the root of the path (i.e. "build", "raw" etc.)

To use this utility you will need to ensure that your AWS credentials
are available via expected environment variables. Please refer to
the AWS boto3 documentation, which describes this: -

    https://boto3.amazonaws.com/v1/documentation/api/latest/guide/quickstart.html

Alan Christie
February 2019
"""

import argparse
import os
import sys

import boto3

# Expected environment variables (that define the bucket)
s3_bucket_env = 'FRAGALYSIS_S3_BUCKET'
s3_archive_bucket = os.environ.get(s3_bucket_env)
if not s3_archive_bucket:
    print('You must define %s', s3_bucket_env)
    sys.exit(1)

parser = argparse.ArgumentParser('Graph Build File Putter')
parser.add_argument('path', metavar='PATH', type=str,
                    help='The path, relative to the "raw" directory'
                         ' in your S3 bucket. e.g. "activity/senp7/build-1"')

args = parser.parse_args()

# The S3 path path...
dst = + args.path + '/'
s3_client = boto3.client('s3')

# List objects at the destination
target = s3_client.list_objects_v2(Bucket=s3_archive_bucket,
                                   Prefix=dst)
if 'KeyCount' in target and target['KeyCount']:
    print('present')
else:
    print('absent')
