#!/bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )";

# Construct matlab scripts path
matlab_scripts_str="/matlab_scripts";
full_matlab_scripts_path="$DIR$matlab_scripts_str";

# Get input arguments
data_dir=$1;
output_directory=$2;
input_ext=$3;
stromule_channel=$4;
chloroplast_channel=$5;

# Do some initial verification
if [ ! -d "$data_dir" ]; then
  echo "Error: data_dir does not exist!";
  echo "  >> $data_dir was entered as the data_dir";
  echo "  >> Usage: ./RunStromuleTracking.sh data_dir output_directory input_ext stromule_channel chloroplast_channel";
  exit
fi

re='^[0-9]+$'
if ! [[ $stromule_channel =~ $re ]] ; then
   echo "Error: stromule_channel is not a number!"; 
   echo "  >> $stromule_channel was entered as the stromule_channel";
   echo "  >> Usage: ./RunStromuleTracking.sh data_dir output_directory input_ext stromule_channel chloroplast_channel";
   exit 
fi

if ! [[ $chloroplast_channel =~ $re ]] ; then
   echo "Error: chloroplast_channel is not a number!"; 
   echo "  >> $chloroplast_channel was entered as the chloroplast_channel";
   echo "  >> Usage: ./RunStromuleTracking.sh data_dir output_directory input_ext stromule_channel chloroplast_channel";
   exit 
fi

if [ ! -d "$full_matlab_scripts_path" ]; then
  echo "Fatal Error: Cannot find MATLAB scripts path!";
  echo "  >> $full_matlab_scripts_path is where this script is looking for MATLAB files";
  echo "  >> Did you move this file to a different location?";
  echo "  >> There should be a directory called matlab_scripts/ in the same directory";
  exit
fi

# Call matlab script with these arguments
echo "Starting MATLAB Stromule Tracking script via SLURM...";
sbatch -N 1 -c 12 --mem=32000 --partition=gpu --account=gpu --wrap="export PATH=$PATH:/opt/MATLAB/R2018b/bin; matlab -nodesktop -nosplash -nodisplay -r \"addpath('$full_matlab_scripts_path'); try, TrackStromules('$data_dir', '$output_directory', '$input_ext', '$stromule_channel', '$chloroplast_channel'); catch e, disp(getReport(e)), exit(7), end, exit\""
exit;