#!/usr/bin/env python3
# coding=utf-8

"""A utility to get the next available build number for a build path.
The output of the utility is the next build number, i.e. '2'.

To use this utility you will need to ensure that your AWS credentials
are available via expected environment variables. Please refer to
the AWS boto3 documentation, which describes this: -

    https://boto3.amazonaws.com/v1/documentation/api/latest/guide/quickstart.html

Alan Christie
March 2019
"""

import argparse
import os
import sys

import boto3

# Expected environment variables (that define the bucket)
s3_bucket_env = 'FRAGALYSIS_S3_BUCKET'
s3_build_root = 'build'
s3_archive_bucket = os.environ.get(s3_bucket_env)
if not s3_archive_bucket:
    print('You must define %s' % s3_bucket_env)
    sys.exit(1)

parser = argparse.ArgumentParser('Graph Build Number Getter')
parser.add_argument('path', metavar='PATH', type=str,
                    help='The path build in your S3 bucket where'
                         ' you expect the build to be placed'
                         ' e.g. "activity/senp7"')

args = parser.parse_args()

s3_client = boto3.client('s3')

# Query S3 until we find a build that's not used.
next_free_number = 1
found = False

# The S3 path path...
dst = s3_build_root + '/' + args.path + '/'

# List objects at the destination.
# This will include things like
# 'build/activity/senp7/build-1/nodes.txt.gz'.
# For each key, look for a 'build' and then inspect its number.
target = s3_client.list_objects_v2(Bucket=s3_archive_bucket,
                                   Prefix=dst)
if 'Contents' in target:
    for content in target['Contents']:
        fields = content['Key'].split('/')
        for field in fields:
            if field.startswith('build-'):
                number = int(field.split('-')[1])
                if number >= next_free_number:
                    next_free_number = number + 1

# The next available build number...
print(next_free_number)
