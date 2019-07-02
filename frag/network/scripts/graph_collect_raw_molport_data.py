#!/usr/bin/env python3
# coding=utf-8

"""A utility to collect raw (vendor) data files from the vendor sites
and upload them to a new raw location on AWS S3. The output of this
program is a string representing the path to any new data on S3,
which will be something like 'vendor/molport/2019-06'. The string is also
written to the file 'collected-release-id.txt' in the collection directory.

If there is no new data the program exit status is non-zero.

There may be output if there's an error, so output must only be interpreted
as a path if the exit code of the program is 0.

This module assumes that the data resides in an S3 bucket with a directory
structure that complies with the **Fragalysis S3 Data Archive** definition.

To use this utility you will need to ensure that your AWS credentials
are available via expected environment variables. Please refer to
the AWS boto3 documentation, which describes this: -

    https://boto3.amazonaws.com/v1/documentation/api/latest/guide/quickstart.html

The collection username is expected to be available in the
MOLPORT_COLLECT_USERNAME and MOLPORT_COLLECT_PASSWORD variables.

Alan Christie
June 2019
"""

import argparse
from ftplib import FTP
import logging
import os
import re
import shutil
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

collect_username = os.environ.get('MOLPORT_COLLECT_USERNAME')
collect_password = os.environ.get('MOLPORT_COLLECT_PASSWORD')
if not collect_username or not collect_password:
    logger.error('You must define a username and password')
    sys.exit(1)

# Molport release directories
# i.e. 2018-07
SET_RE = re.compile(r'^\d\d\d\d-\d\d$')
# The FTP site
FTP_HOME = 'molport.com'
FTP_ROOT = 'molport_database'
S3_STORAGE_PATH = 'vendor/molport'

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
    """An FTP callback method used by check_latest().
    Used to find the latest release on the customer site.
    If a release is found, the most recent is sets int the global
    variables latest_release_str (i.e. the string '2019-06') and
    latest_release_id (i.e. the number 201906)
    """
    global latest_release_id
    global latest_release_str

    # A MolPort file looks something
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


def ftp_record_callback(context):
    """An FTP callback method used by get_latest().
    Used to put the release files (the SMILES files)
    into a global list (latest_release_files) which is then
    processed from the get_latest() method.
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


def get_latest():
    """Collects all the files for the release identified by the
    'latest_release_str' string's value.
    """
    if not latest_release_str:
        return

    ftp = FTP(FTP_HOME)
    ftp.login(collect_username, collect_password)

    # The release files we want are in the FTP ROOT
    # at '<release>/All Stock Compounds/SMILES'.
    ftp.cwd(os.path.join(FTP_ROOT,
                         latest_release_str,
                         'All Stock Compounds',
                         'SMILES'))
    ftp.retrlines('LIST', ftp_record_callback)

    # Now we've listed the directory,
    # using the current FTP instance,
    # retrieve the files...
    for filename in latest_release_files:
        dest_path = os.path.join(collect_dir, filename)
        ftp.retrbinary('RETR %s' % filename, open(dest_path, 'wb').write)

    ftp.quit()

    # Upload the list of files (removing the source as we go)...
    s3_client = boto3.client('s3')
    for filename in latest_release_files:
        src = os.path.join(collect_dir, filename)
        dst = s3_data_root + '/' + S3_STORAGE_PATH + '/' + latest_release_str + '/' + filename
        s3_client.upload_file(src, s3_archive_bucket, dst)
        os.remove(src)


def check_held():
    """Gets the latest data held on S3, setting latest_held_id
    and latest_held_str.
    """
    global latest_held_id
    global latest_held_str

    src = s3_data_root + '/' + S3_STORAGE_PATH + '/'
    # Get the contents of the selected directory
    # and look for potential releases...
    s3_client = boto3.client('s3')
    resp = s3_client.list_objects_v2(Bucket=s3_archive_bucket,
                                     Prefix=src)
    for content in resp['Contents']:
        content_parts = content['Key'].split('/')
        if len(content_parts) >= 4:
            potential_collection = content_parts[3]
            if SET_RE.match(potential_collection):
                held_id = to_release_id(potential_collection)
                if held_id > latest_held_id:
                    latest_held_id = held_id
                    latest_held_str = potential_collection


parser = argparse.ArgumentParser('Graph MolPort File Collector')
parser.add_argument('collection', metavar='DIR', type=str,
                    help='The local collection directory for any'
                         ' collected data, which must not exist and'
                         ' will be created')
parser.add_argument("--force", dest="force", action="store_true")

args = parser.parse_args()
collect_dir = args.collection

# Destination?
if not args.force and os.path.isdir(collect_dir):
    logger.error('Collection exists')
    sys.exit(1)

# Wipe any existing collection directory
if os.path.isdir(collect_dir):
    shutil.rmtree(collect_dir)
# Create the destination directory
if not os.path.isdir(collect_dir):
    os.makedirs(collect_dir)

# Find the most recent collected data.

check_held()

# Now inspect the vendor site to see the most recent data.
#
# If we've missed some we just want the very latest, not the next logical.
# So if we have data for June and new data exists for July and August
# we want August's data.

check_latest()

# If there's no new data: -
#
#   exit with a non-zero code
#
# Otherwise: -
#
# - Download the new data to the collection directory
# - Upload to the S3 path
# - Print the location of the new data.

if latest_release_id <= latest_held_id:
    print('Nothing newer than held release found')
    print('Latest held is %s (%s)' % (latest_held_str, latest_held_id))
    print('Latest release is %s (%s)' % (latest_release_str, latest_release_id))
    sys.exit(2)

# There's new data.
# Download and store on S3
get_latest()

# Finally, print the new path,
# this is the 'path' argument with the new raw directory appended,
# and will be something like 'vendor/molport/2019-06'
collected_release_id = S3_STORAGE_PATH + '/' + latest_release_str
print(collected_release_id)

# And write the collected release ID to file
# (to simplify inter-playbook execution)
id_file = open(os.path.join(collect_dir, 'collected-release-id.txt'), 'wt')
id_file.write(collected_release_id)
id_file.close()
