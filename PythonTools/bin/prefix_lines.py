#!/common/software/install/migrated/anaconda/python3-2020.07-mamba/bin/python

"""
Simple script to add a prefix to every line in a file.
"""

import sys

for line in sys.stdin:
    print(sys.argv[1] + line, end=" ")
