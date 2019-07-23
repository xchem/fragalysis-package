#!/usr/bin/env python3
# coding=utf-8

"""A utility to get the next available combination number.
The output of the utility is the next combination number, i.e. '2'.

To use this utility you will need to ensure that your AWS credentials
are available via expected environment variables. Please refer to
the AWS boto3 documentation, which describes this: -

    https://boto3.amazonaws.com/v1/documentation/api/latest/guide/quickstart.html

Alan Christie
March 2019
"""

import os
import sys

import boto3

# Expected environment variables (that define the bucket)
s3_bucket_env = 'AWS_S3_BUCKET'
s3_combination_root = 'combination'
s3_archive_bucket = os.environ.get(s3_bucket_env)
if not s3_archive_bucket:
    print('You must define %s' % s3_bucket_env)
    sys.exit(1)

s3_client = boto3.client('s3')

# Query S3 until we find a build that's not used.
next_free_number = 1

# The S3 path path...
dst = s3_combination_root + '/'

# List objects at the destination.
# This will include things like
# 'combination/1'.
target = s3_client.list_objects_v2(Bucket=s3_archive_bucket,
                                   Prefix=dst)
if 'Contents' in target:
    for content in target['Contents']:
        fields = content['Key'].split('/')
        for field in fields:
            is_combination_number = True
            try:
                number = int(field)
            except ValueError:
                is_combination_number = False
            if is_combination_number:
                if number >= next_free_number:
                    next_free_number = number + 1

# The next available combination number...
print(next_free_number)
