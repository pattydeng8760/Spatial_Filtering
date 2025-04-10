# Function to run the filtering process
from .spatial_filter import SpatialFilter
import argparse

# The required inputs
dir_to_post='/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/Isosurface/'          # The directory to the source file
mesh_Dir = '/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/MESH_ZONE_Apr24/'               # The mesh directory
mesh_fileName = 'Bombardier_10AOA_Combine_Apr24.mesh.h5'                                                # The source mesh file name in .h5 format
cut_location = 'Cut_5mm_tip_VGT'                                                            # The source file name
variable = ['dilatation']                                                           # The variable of interest


# Define test parameters
args = argparse.Namespace(
    dir_to_post=dir_to_post,
    mesh_dir=mesh_Dir,
    mesh_filename=mesh_fileName,
    cut_location=cut_location,
    dt=0.155e-7 * 2000,
    freq_min=800,
    freq_max=1500,
    zones=5000,
    variable=variable,
    cores=20,
    override=False,
    nstart = 200,
    nend = 600,
    ndt = 2,
)

# Call the FieldFilter class like you would from the CLI
filter_job = SpatialFilter(args)
filter_job.run()