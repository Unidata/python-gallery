#!/usr/bin/env python

import fileinput
from pathlib import Path
import sys

# Recursively grab examples
scripts = Path('examples').rglob('*.py')

# Go through each script, looking for http://thredds.ucar.edu and replacing
# it with http://thredds-test.unidata.ucar.edu or
# http://thredds-test.unidata.ucar.edu to test against

search_text = 'http://thredds.ucar.edu'

if sys.argv[1] == 'main':
    sys.exit(0)
elif sys.argv[1] == 'test':
    replace_text = 'http://thredds-test.unidata.ucar.edu'
elif sys.argv[1] == 'dev':
    replace_text = 'http://thredds-dev.unidata.ucar.edu'
elif sys.argv[1] == 'atm':
    replace_text = 'http://atm.ucar.edu'
elif sys.argv[1] == 'jetstream':
    replace_text = 'http://thredds-jetstream.unidata.ucar.edu'
else:
    raise ValueError('Invalid THREDDS specifier: {}'.format(sys.argv[1]))

print('Converting examples to use {}'.format(replace_text))

for script in scripts:
    print('Processing {}'.format(str(script)))
    with fileinput.FileInput(str(script), inplace=True) as f:
        for line in f:
            print(line.replace(search_text, replace_text), end='')
