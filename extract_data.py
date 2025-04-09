# data_extraction.py
import h5py
import os
import numpy as np
import glob

def extract_data(dir_to_post:str, mesh:str, output:str, var_interest:str, reload_data:bool=False):
    """Extract the cutplane data from AVBP cutplanes to a single .h5 file for each variable of interest.
    Args:
        dir_to_post (str): The directory where the AVBP cutplane files are located
        mesh (str): The mesh file path for the cutplane, this is extracted though the mesh_utils.py
        output (str): The destination directory where the extracted data is to be stored
        var_interest (list): The list of variables to be extracted from the AVBP cutplane files
        reload_data (bool, optional): If True, reloads the data even if it already exists. Defaults to False.
    Returns:
        nnodes (int): The number of nodes in the data
        nb_times (int): The number of timesteps extracted
    """
    # Print header for base data extraction
    header = "Extracting Base Data"
    print(f"\n{header:.^80}\n")
    
    # Get a sorted list of all .h5 files in the given directory
    file_list = sorted(glob.glob(os.path.join(dir_to_post, '*.h5')))
    nb_times = len(file_list)
    print(f"The number of timesteps extracted is: {nb_times}.")
    
    # Open the first file using a context manager to determine the number of nodes
    with h5py.File(file_list[0], 'r') as fid:
        # Using a fixed key path to get the node count
        nnodes = np.shape(fid['0000/instants/0000/variables/dilatation_node'])[0]
    print(f"The number of nodes in the data is: {nnodes}.")
    
    # Loop over each variable of interest
    for var in var_interest:
        var_name = var.split('/')[-1]
        if os.path.exists(os.path.join(output, f"field_{var_name}.h5")) and reload_data == False:
            print(f"File {var_name} already exists. Skipping extraction.")
            return nnodes, nb_times
        else:
            print(f"\nExtracting var: {var}")
            # Initialize an array to hold the data for this variable across all timesteps
            data = np.zeros((nnodes, nb_times))
            it = 0  # Initialize our counter for the while loop
            # Read each file with a while loop, ensuring file closure with a context manager
            while it < nb_times:
                filename = file_list[it]
                with h5py.File(filename, 'r') as fid:
                    data[:, it] = fid[var][()]
                # Provide progress feedback every 100 files and on the last file
                if it % 100 == 0 or it == nb_times - 1:
                    print(f"    {it} | {nb_times - 1} files loaded")
                it += 1
            # Save the collected data for the current variable
            field_file = os.path.join(output, f"field_{var_name}.h5")
            with h5py.File(field_file, 'w') as gid:
                gid.create_dataset('field', data=data)
            print(f"\nThe base data is saved in: {field_file}")
    # Print footer to indicate the extraction is complete
    footer = "Base Data Extraction Complete"
    print(f"\n{footer:.^80}\n")
    return nnodes, nb_times