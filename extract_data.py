# data_extraction.py
from antares import *
import h5py
import os
import numpy as np
import glob

def extract_data(dir_to_post:str, cut_location:str ,mesh:str, output:str, var_interest:str, data_type:str, reload_data:bool=False):
    """Extract the cutplane data from AVBP cutplanes to a single .h5 file for each variable of interest.
    Args:
        dir_to_post (str): The directory where the AVBP cutplane files are located
        cut_location (str): The location of the cut plane, defined either explicitly via PIV planes
        mesh (str): The mesh file path for the cutplane, this is extracted though the mesh_utils.py
        output (str): The destination directory where the extracted data is to be stored
        var_interest (list): The list of variables to be extracted from the AVBP cutplane files
        data_type (str): The type of data to be extracted, options are 'FWH', 'OUTPUT', 'EXTRACT', or 'CLIP'
        reload_data (bool, optional): If True, reloads the data even if it already exists. Defaults to False.
    Returns:
        nnodes (int): The number of nodes in the data
        nb_times (int): The number of timesteps extracted
    """
    # Print header for base data extraction
    header = "Extracting Base Data"
    print(f"\n{header:.^60}\n")
    assert data_type in ['FWH', 'OUTPUT', 'EXTRACT', 'CLIP'], "Invalid data type. Choose from 'FWH', 'OUTPUT', 'EXTRACT', or 'CLIP'."
    # Get a sorted list of all .h5 files in the given directory
    file_list = sorted(glob.glob(os.path.join(dir_to_post,cut_location, '*.h5')))
    nb_times = len(file_list)
    print(f"The number of timesteps extracted is: {nb_times}.")
    
    # Open the first file using a context manager to determine the number of nodes
    with h5py.File(file_list[0], 'r') as fid:
        # Using a fixed key path to get the node count
        nnodes = np.shape(fid['0000/instants/0000/variables/dilatation_node'])[0]
    print(f"The number of nodes in the data is: {nnodes}.")
    
    # Loop over each variable of interest
    for var in var_interest:
        if os.path.exists(os.path.join(output, f"field_{var}.h5")) and reload_data == False:
            print(f"File {var} already exists. Skipping extraction.")
            header = "Base Data Extracted"
            print(f"\n{header:.^60}\n")
            return nnodes, nb_times
        else:
            print(f"\n---->Extracting var: {var}")
            # Initialize an array to hold the data for this variable across all timesteps
            data = np.zeros((nnodes, nb_times))
            for it, filename in enumerate(file_list):
                if data_type == 'FWH':
                    f=h5py.File(filename,'r')
                    press=f['frame_data/'+var][()]
                    data[:,it]=press
                    f.close()
                elif data_type == 'OUTPUT':
                    r = Reader('hdf_avbp')
                    r['filename'] = filename
                    base  = r.read()
                    press = base[0][0][(var, 'node')]
                    data[:,it]=press
                elif data_type == 'EXTRACT' or data_type == 'CLIP':
                    r = Reader('hdf_antares')
                    r['filename'] = filename
                    base  = r.read()
                    press = base[0][0][(var, 'node')]
                    data[:,it]=press
                # Provide progress feedback every 100 files and on the last file
                if it % 100 == 0 or it == nb_times - 1:
                    print('    {0:3d} | {1:3d} files loaded'.format(it, nb_times-1))
                    
            # Save the collected data for the current variable
            field_file = os.path.join(output, f"field_{var}.h5")
            with h5py.File(field_file, 'w') as gid:
                gid.create_dataset('field', data=data)
            print(f"\n---->The base data is saved in: {field_file}")
            
    # Print footer to indicate the extraction is complete
    footer = "Base Data Extraction Complete"
    print(f"\n{footer:.^60}\n")
    return nnodes, nb_times