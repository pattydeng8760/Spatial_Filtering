#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=0-06:30:00
#SBATCH --job-name=Filtering
#SBATCH --mail-user=patrickgc.deng@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-plavoie

source /project/m/moreaust/Env/avbpNpy_env.sh
use_py_tools
export PYTHONPATH="/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/PostProc/Filtering:$PYTHONPATH"

# Define the static paths and parameters
DIR_TO_POST="/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/Isosurface/"
MESH_DIR="/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/MESH_ZONE_Apr24/"
MESH_FILENAME="Bombardier_10AOA_Combine_Apr24.mesh.h5"
# dt is computed as: 0.155e-7 * 2000
DT=$(python -c "print(0.155e-7 * 2000)")

# Other fixed parameters
ZONES=5000
VARIABLE="dilatation"
CORES=20
OVERRIDE="False"
NSTART=200
NEND=600
NDT=2

# Define arrays for the planes (cut locations)
planes=("Cut_5mm_tip_VGT" "Cut_5mm_tip_VGT" "Cut_midspan_TE_VGT")
# Define arrays for frequency ranges; here, each frequency range is defined by its min and max values.
# Make sure both arrays have the same length.
freq_mins=(100 800 1500 2500)
freq_maxs=(1000 2000 2500 4500)

# Loop over each plane
for plane in "${planes[@]}"; do
    # Loop over the indices of the frequency arrays
    for i in "${!freq_mins[@]}"; do
        freq_min=${freq_mins[i]}
        freq_max=${freq_maxs[i]}
        
        echo "Running filter for plane: ${plane} with frequency range: ${freq_min} - ${freq_max}"
        
        # Run the Python script with the desired parameters
        python -m filtering.spatial_filter \
            --dir_to_post "${DIR_TO_POST}" \
            --mesh_dir "${MESH_DIR}" \
            --mesh_filename "${MESH_FILENAME}" \
            --cut_location "${plane}" \
            --dt "${DT}" \
            --freq_min "${freq_min}" \
            --freq_max "${freq_max}" \
            --zones "${ZONES}" \
            --variable "${VARIABLE}" \
            --cores "${CORES}" \
            --override "${OVERRIDE}" \
            --nstart "${NSTART}" \
            --nend "${NEND}" \
            --ndt "${NDT}"
            
        # Optionally add a sleep interval to avoid overlapping jobs
        # sleep 1
    done
done
