#!/bin/bash
#PBS -l select=1:ncpus=16:mem=128gb
#PBS -l walltime=72:00:00
#PBS -N ST_pipeline
#PBS -j oe

set -euo pipefail

cd "$PBS_O_WORKDIR"

#########################################
##                                     ##
## Output some useful job information. ##
##                                     ##
#########################################

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: current home directory is $PBS_O_HOME
echo ------------------------------------------------------

echo "Starting setup..."

# Load Conda via Miniforge 
echo "Initializing Conda..."
eval "$(~/miniforge3/bin/conda shell.bash hook)"

# Check if the environment exists; if not, create it
if ! conda env list | grep -q "ST_env"; then
    echo "Environment ST_env not found. Creating it from environment.yml..."
    conda env create -f environment.yml
else
    echo "ST_env already exists. Skipping creation..."
fi

# Activate the environment 
echo "Activating ST_env..."
conda activate ST_env

# install specific ver of setuptools
echo "Installing setuptools..."
pip install --upgrade setuptools==81.0.0 wheel

# Install muspan 
MUSPAN_USER="GetMuSpAn"
MUSPAN_PASS="SpatialBiology"

echo "Downloading and installing muspan..."
curl -L -u "${MUSPAN_USER}:${MUSPAN_PASS}" -o muspan.zip "https://docs.muspan.co.uk/code/latest.zip"
pip install muspan.zip
rm muspan.zip

# Install requirements
echo "Installing requirements.txt..."
pip install --no-cache-dir -r requirements.txt

# Run pipeline
echo "Environment is ready! Running pipeline..."
python Modules/__main__.py config_Xenium.toml 

echo "Pipeline finished successfully!"

# Deactivate
conda deactivate