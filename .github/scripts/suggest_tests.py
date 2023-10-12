#!/usr/bin/env python

import os
import argparse
import glob
import re

rootdir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
testdir = os.path.join(rootdir, 'test')

parser = argparse.ArgumentParser()
parser.add_argument("file", nargs='+', help="list of files that were changed")
args = parser.parse_args()

# find all test scripts
tests = glob.glob(os.path.join(testdir, '*.m'))

# this will hold a list of test files with either public or private data
data_public = []
data_private = []

# this will hold the test scripts that have a dependency on any of the inputs
suggestions = []

for testfile in tests:
    pattern1 = '^% DATA.*\\bpublic|no\\b'  # this includes NO DATA
    pattern2 = '^% DATA.*\\bprivate\\b'

    for line in open(os.path.join(testdir, testfile)):
        if re.search(pattern1, line):
            data_public.append(testfile)
        if re.search(pattern2, line):
            data_private.append(testfile)

for changedfile in args.file:
    f, x = os.path.splitext(os.path.basename(changedfile))
    pattern = '^% DEPENDENCY.*\\b' + f + '\\b'

    for testfile in tests:
        for line in open(os.path.join(testdir, testfile)):
            if re.search(pattern, line):
                suggestions.append(testfile)

# remove duplicates
data_public = set(data_public)
data_private = set(data_private)
suggestions = set(suggestions)

suggestions_public = suggestions.intersection(data_public)
suggestions_private = suggestions.intersection(data_private)

if len(suggestions_public) or len(suggestions_private):
    print('You should test whether your modifications do not break anything.')
    print('See <https://www.fieldtriptoolbox.org/development/testing/>')
    print()

if len(suggestions_public):
    print('When outside the DCCN, please consider testing: ', end='')
    for i, file in enumerate(suggestions_public):
        f, x = os.path.splitext(os.path.basename(file))
        if i==len(suggestions_public)-1:
            print(f, end='\n')
        else:
            print(f, end=', ')
    print()

if len(suggestions_private):
    print('When inside the DCCN, please also consider testing: ', end='')
    for i, file in enumerate(suggestions_private):
        f, x = os.path.splitext(os.path.basename(file))
        if i==len(suggestions_private)-1:
            print(f, end='\n')
        else:
            print(f, end=', ')
    print()

if len(suggestions_public) and len(suggestions_private):
    print('Suggested tests outside the DCCN use public data or do not use data.')
    print('Suggested tests inside the DCCN use private data.') 