import os
import json
import shutil
import subprocess
import numpy as np
from itertools import product

import glob
import re

# Get matching systems
pattern = "system_*.json"
system_list = glob.glob(pattern)
print("Files:", system_list)
temp_list = []
pattern = re.compile(r"system_(\d+\.\d+)\.json")

for filename in system_list:
    match = pattern.match(filename)
    temp = match.group(1)

# Create a base directory for the parameter sweep results
base_dir = '$SCRATCH/metathesis_sweep_01'
os.makedirs(base_dir, exist_ok=True)
# Construct the path to the executable (passed to sbatch, main is large, so we don't want thousands of copies)
os.environ['main_executable'] = os.path.join(os.path.abspath(os.curdir), 'main')

# Loop over each system file
for filename in system_list:
    match = pattern.match(filename)
    temp = match.group(1)
    # Create a directory for this parameter combination
    unique_name = f'T={temp}'
    dir_name = f'{base_dir}/{unique_name}'
    os.makedirs(dir_name, exist_ok=True)

    # Copy other necessary files to the new directory
    shutil.copy('submit_sweep.sh', dir_name)
    shutil.copy('parameters.prm', dir_name)
    shutil.copy('BC_AMR.prm', dir_name)
    shutil.copy('plot_and_save.py', dir_name)
    shutil.copy(filename, f'{dir_name}/system.json')

    # Run 'main' inside the new directory
    subprocess.run(['sbatch', \
                    '-J', unique_name, \
                    '-o', f'{unique_name}.out', \
                    '-e', f'{unique_name}.err', \
                    'submit_sweep.sh'], \
                    cwd=dir_name, env=os.environ)

print('Parameter sweep submission completed.')
