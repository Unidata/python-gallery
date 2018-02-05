#!/usr/bin/env python

import fileinput
from pathlib import Path


# Recursively grab examples
scripts = Path('examples').rglob('*.py')

# Go through each script, looking for http://thredds.ucar.edu and replacing
# it with http://thredds-test.unidata.ucar.edu to test against
search_text = 'http://thredds.ucar.edu'
replace_text = 'http://thredds-test.unidata.ucar.edu'

for script in scripts:
    print('Processing {}'.format(str(script)))
    with fileinput.FileInput(str(script), inplace=True, backup='.bak') as f:
        for line in f:
            print(line.replace(search_text, replace_text), end='')
