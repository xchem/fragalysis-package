#!/usr/bin/env python3
# coding=utf-8

"""A utility to collect raw (vendor) data files from the vendor sites
and upload them to a new raw location on AWS S3. The output of this
program is a string representing the path to any new data omn S3. If there is
no new data the program output is blank. There will be output if
there's an error, so output must only be interpreted as a path if the
exit code of the program is 0.

This module assumes that the data resides in an S3 bucket with a directory
structure that complies with the **Fragalysis S3 Data Archive** definition.

To use this utility you will need to ensure that your AWS credentials
are available via expected environment variables. Please refer to
the AWS boto3 documentation, which describes this: -

    https://boto3.amazonaws.com/v1/documentation/api/latest/guide/quickstart.html

The collection username is expected to be available in the
FRAGALYSIS_COLLECT_USERNAME and FRAGALYSIS_COLLECT_PASSWORD variables.

Alan Christie
June 2019
"""

import argparse
from ftplib import FTP
import logging
import os
import re
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
s3_data_root = 'raw'
s3_archive_bucket = os.environ.get(s3_bucket_env)
if not s3_archive_bucket:
    logger.error('You must define %s', s3_bucket_env)
    sys.exit(1)

collect_username = os.envrion.get('FRAGALYSIS_COLLECT_USERNAME')
collect_password = os.envrion.get('FRAGALYSIS_COLLECT_PASSWORD')
if not collect_username or not collect_password:
    logger.error('You must define a username and password')
    sys.exit(1)

# Molport release directories
# i.e. 2018-07
SET_RE = re.compile(r'^\d\d\d\d-\d\d$')
# The FTP site
FTP_HOME = 'molport.com'
FTP_ROOT = 'molport_database'
S3_STORAGE_PATH = 'raw/vendor/molport'

latest_release_str = None
latest_release_id = 0
latest_release_files = []
latest_held_str = None
latest_held_id = 0
collect_dir = None


def to_release_id(release):
    """Turns a release string into an ID.
    Call with a release like '2017-01'
    """
    return int(release[:4] + release[5:7])


def ftp_sniff_callback(context):
    """An FTP callback method.
    """
    global latest_release_id
    global latest_release_str

    # A MolPort file looks something
    # like this as a directory listing: -
    #
    # drwxr-xr-x   6 root     root         4096 Aug  6  2018 2018-07
    parts = context.split()
    if len(parts) != 9 or parts[0][0] != 'd' or not SET_RE.match(parts[8]):
        return
    # Convert to a numeric identity,
    # high values are later releases.
    rel_id = to_release_id(parts[8])
    if rel_id > latest_release_id:
        latest_release_id = rel_id
        latest_release_str = parts[8]


def ftp_pull_callback(context):
    """An FTP callback method.
    """
    global latest_release_id
    global latest_release_str

    # A MolPort file looks something
    # like this as a directory listing: -
    #
    # drwxr-xr-x   6 root     root         4096 Aug  6  2018 2018-07
    parts = context.split()
    if len(parts) == 0 or not parts[-1].endswith('gz'):
        return
    # Convert to a numeric identity,
    # high values are later releases.
    filename = parts[-1]
    latest_release_files.append(filename)


def check_latest():
    """Iterates through the release directory.
    After this call 'latest_release_id' is set according to the latest release
    and 'latest_release_str' is the release directory.
    """
    ftp = FTP(FTP_HOME)
    ftp.login(collect_username, collect_password)
    ftp.cwd(FTP_ROOT)
    ftp.retrlines('LIST', ftp_sniff_callback)
    ftp.quit()

    logger.info('Latest release is %s', latest_release_str)

def collect():
    """Collects all the files for the release identified by the
    'latest_release_str' string's value.
    """
    if not latest_release_str:
        return

    ftp = FTP(FTP_HOME)
    ftp.login(collect_username, collect_password)
    ftp.cwd(os.path.join(FTP_ROOT, 'All Stock Compounds', 'SMILES'))
    ftp.retrlines('LIST', ftp_pull_callback)

    # Now retrieve the files...
    for filename in latest_release_files:
        logger.info('Getting %s...', filename)
        ftp.retrbinary('RETR %s' % filename, open(filename, 'wb').write)

    ftp.quit()

    # Upload the list of files...
#    s3_client = boto3.client('s3')
    potential_files = os.listdir(collect_dir)
    for potential_file in potential_files:
        src = os.path.join(args.source, potential_file)
        if os.path.isfile(src):
            dst = S3_STORAGE_PATH + '/' + latest_release_str + '/' + potential_file
            logger.info('Putting %s -> %s...', potential_file, latest_release_str)
#            s3_client.upload_file(src, s3_archive_bucket, dst)


def check_held():
    """Gets the latest held data, setting latest_held_id.
    """
    global latest_held_id
    global latest_held_str

    src = S3_STORAGE_PATH + '/'
    # Get the contents of the selected directory
    # and look for potential releases...
    s3_client = boto3.client('s3')
    resp = s3_client.list_objects_v2(Bucket=s3_archive_bucket,
                                     Prefix=src)
    for content in resp['Contents']:
        if content['Key'].endswith('/'):
            potential_collection = content['Key'].split('/')[-2]
            if SET_RE.match(potential_collection):
                held_id = to_release_id(potential_collection)
                if held_id > latest_held_id:
                    latest_held_id = held_id
                    latest_held_str = potential_collection

    logger.info('Latest held is %s', latest_held_str)


parser = argparse.ArgumentParser('Graph MolPort File Collector')
parser.add_argument('collection', metavar='DIR', type=str,
                    help='The local collection directory for any'
                         ' collected data, which must not exist and'
                         ' will be created')

args = parser.parse_args()
collect_dir = args.collection

# Destination?
if not args.force and os.path.isdir(collect_dir):
    logger.error('Collection exists')
    sys.exit(1)

# Find the most recent collected data.

check_held()

# Now inspect the vendor site to see the most recent data.
#
# If we've missed some we just want the very latest, not the next logical.
# So if we have data for June and new data exists for July and August
# we want August's data.

check_latest()

# If there's new data: -
#
# - Download it to the collection directory
# - Upload to the S3 path

if latest_release_id > latest_held_id:

    # There's new data.
    # Download and store on S3
    collect()

    # Finally, print the new path,
    # this is the 'path' argument with the new raw directory appended
    print(S3_STORAGE_PATH + '/' + latest_release_str)
