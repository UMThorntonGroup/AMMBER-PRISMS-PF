#!/bin/bash

#SBATCH -A DMR110007
#SBATCH -p wholenode
#SBATCH -N 2
#SBATCH --ntasks=256
#SBATCH -t 01:00:00
#SBATCH --mail-user=xmen@umich.edu
#SBATCH --mail-type=FAIL

# run this using sbatch submit_sweep.sh 
# run from the results directory or with cwd set to the results directory

module load ffmpeg
#module load visit
pwd
date

# Run simulation
mpirun -np $SLURM_NTASKS $main_executable
# Generate frames and movie
visit -cli -nowin -s plot_and_save.py
ffmpeg -framerate 30 -i frames/frame_%04d.png -c:v mpeg4 -pix_fmt yuv420p -q:v 2 movie_$SLURM_JOB_NAME.mp4

# Remove excess output files
# Loop through all files starting with 'solution-' and ending with '.vtu'
# Sort them numerically (assuming filenames have a numeric part)
#count=0
#find . -maxdepth 1 -type f -regex './solution-.*\.vtu' | sort | while read -r file; do
#    # Delete the file if it's not every 100th file
#    if ((count % 100 != 0)); then
#        rm "$file"
#    fi
#    # Increment counter
#    ((count++))
#done

# Make a zip of everything for safekeeping / quick file transfer
zip ../$SLURM_JOB_NAME.zip *
zip ../movies.zip movie_$SLURM_JOB_NAME.mp4
last_frame=$(ls -v frames/frame_*.png | tail -n 1)
cp $last_frame ../final_$SLURM_JOB_NAME.png
zip -j ../final_frames.zip ../final_$SLURM_JOB_NAME.png
