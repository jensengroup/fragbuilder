## FragBuild 0.1 Copyright (c) 2012 Anders S. Christensen
## Examples
import fragbuild
import string
import math
import sys

log_filename = sys.argv[1]

# Load results
nmr_data = fragbuild.nmr_file(log_filename)

# Print all CA shielding values
print nmr_data.get_chem_shifts("XX")

