#!/usr/bin/env python3

'''
C. Lechner, European XFEL, 2024-01-05

Function to obtain ID of most recent git commit.

Note that this function does not report any uncommitted changes.

Example:
>>> from ocelot.utils.git_commit_id import get_git_commit_id
>>> x=get_git_commit_id()
>>> print(x)
16425bd808d2c45e62d512093ee51190969fa5d7
>>>
'''

import subprocess
import os

def get_git_commit_id():
    # To get a directory within the active OCELOT directory hierarchy,
    # extract directory name from filename of this module
    wd=os.path.dirname(__file__)
    try:
        cp = subprocess.run(['git', 'log', '--format=%H', '-n',  '1'], 
            capture_output=True, cwd=wd)
    except FileNotFoundError:
        print('Error: Unable to obtain git commit id')
        return None
    
    if cp.returncode:
        print('Error: Running git gave non-zero exit code')
        return None

    # Command is expected to return a single line terminated by '\n': the id of the most recent commit
    commitid = cp.stdout
    commitid = commitid.decode().splitlines()
    commitid = commitid[0]
    return(commitid)

