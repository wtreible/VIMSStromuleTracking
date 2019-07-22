#!/usr/bin/env python
from __future__ import print_function
import os, sys, argparse, shlex, subprocess

verbose = True

#================================#
# ==== User Input Arguments ==== #
#================================#
def get_args():
  parser = argparse.ArgumentParser(description='VIMS Lab Stromule Tracking v0.10')

  parser.add_argument("--input-path", type=str, dest="input_path", required=True,
                      help="Path to the input timeseries segmentation")   
  parser.add_argument("--output-path", type=str, dest='output_path', required=True,
                      help="Path to the output directory.")
  parser.add_argument("--input-ext", type=str, dest='input_ext', default='.tif',
                      help="Extension of the segmentation images (Default: .tif)")             
  parser.add_argument("--stromule-channel", type=int, dest='stromule_channel', default=1,
                      help="Image channel that stromules are on in the segmented image (Default: 1; Channels start at 1).")
  parser.add_argument("--chloroplast-channel", type=int, dest='chloroplast_channel', default=2,
                      help="Image channel that chloroplasts are on in the segmented image (Default: 2; Channels start at 1).")
  return parser.parse_args()

  
# Parse Params
args = get_args()
input_path = args.input_path
output_path = args.output_path
input_ext = args.input_ext
stromule_channel = args.stromule_channel
chloroplast_channel = args.chloroplast_channel

#==================================#
# ==== Verification of Inputs ==== #
#==================================#

# Script paths
self_path = os.path.realpath(__file__)
self_dir = os.path.dirname(self_path)
matlab_scripts_path = os.path.join(self_dir, 'matlab_scripts') 

# Verify all paths exist before submitting slurm job
if not os.path.isdir(matlab_scripts_path):
  sys.exit("Error: The matlab_scripts directory could not be found (this script should be in the same folder as matlab_scripts)!")

if not (os.path.isdir(input_path)):
  sys.exit("Error: input_path directory does not exist!")
  
if not os.path.isdir(output_path):
  print("Warning: output_path directory does not exist; this directory will be created... (%s)" % (output_path,))
  # Verify the containing directory exists for the folder to be made...
  # (i.e. if /a/b/c doesn't exist, verify /a/b exists to make dir c)
  up_one_output_path = os.path.dirname(output_path if output_path[-1] != '/' else output_path[:-1])
  if not os.path.isdir(up_one_output_path):
    sys.exit("Error: The directory containing output_path (%s) doesn't exist either! Exiting..." % (up_one_output_path,))

#===============================#
# ==== Build Slurm Command ==== #
#===============================#

# Define debug print function if verbose mode is on
debug_print = print if verbose else lambda *args: None

# Initialize Sbatch Command
sbatch_flags = ["-N 1",
                "-c 24",
                "--mem=32000",
                "--partition=gpu",
                "--account=gpu"]
  
debug_print ("============================")
debug_print ("Sbatch Flags:")
debug_print ("\n".join(sbatch_flags))    

# Construct Matlab Command
matlab_command = "export PATH=$PATH:/opt/MATLAB/R2018b/bin; "
matlab_command += "matlab -nodesktop -nosplash -nodisplay -r \\\"addpath('" + matlab_scripts_path + "'); try, "
matlab_command += "TrackStromules('" + input_path + "', '" + output_path+ "', '" + input_ext + "', '" + str(stromule_channel) + "', '" + str(chloroplast_channel) + "'); "
matlab_command += "catch e, disp(getReport(e)), exit(7), end, exit\\\""
matlab_command_wrapped = '--wrap="' + matlab_command + '"'

# Construct sbatch command
sbatch_command = ["sbatch"] + sbatch_flags + [matlab_command_wrapped]
sbatch_command_str = " ".join(sbatch_command)
debug_print ("============================")
debug_print ("Full SBatch command:")
debug_print (sbatch_command_str)
debug_print ("============================")

final_sbatch_command = shlex.split(sbatch_command_str)
p = subprocess.Popen(final_sbatch_command)
p.wait()
debug_print ("Success!")
debug_print ("============================")
sys.stdout.flush()
