from antares import *
import numpy as np
import os
import sys
import builtins
import shutil
import h5py

# The data reconstrction script
# The data reconstrction script
def reconstruction(mesh,var_interest,nb_times,freq_min:int,freq_max:int,fft_file:str,target_dir:str,nstart=None,nend=None,ndt=None):
    """Applying the reconstruction process to the filtered data stored in fft_file after fft. Loads the mesh file and 
    reconstructs the filtered data at each time step. The reconstructed data is saved in the target_dir.
    Args:
        mesh (string): The mesh file path to be used for the reconstruction
        var_interest (list): The list of variables to be reconstructed
        nb_times (int): The number of time steps in the data
        freq_min (int): The minimum frequency of the filtered data
        freq_max (int): The maximum frequency of the filtered data
        fft_file (string): The path to the filtered data file
        target_dir (string): The directory where the reconstructed data will be saved
        nstart (int): The starting time step for the reconstruction
        nend (int): The ending time step for the reconstruction
        ndt (int): The time step interval for the reconstruction
    Returns:
        None
    """
    # Removing the directory if exists
    if os.path.exists(target_dir):
        shutil.rmtree(target_dir)
    os.makedirs(target_dir, exist_ok = False)        # The mesh file of the transient signal
    text = 'Beginning reconstruction process'
    print(f'\n{text:.^80}\n')  
    # Reading the mesh file
    r=Reader('hdf_antares')
    r['filename']='{0:s}'.format(mesh)
    plane=r.read()
    timearray=[str(i).zfill(4) for i in range(int(nb_times))]
    timearray=timearray[nstart:nend:ndt] if ndt is not None and nend is not None and nstart is not None else timearray
    print("Reading the Base Mesh File")
    print("Reconstructing the filtered signal")
    for var in var_interest:
        # Loading the filtered data from the fft_file
        f=h5py.File(fft_file,'r')
        tim=f['/time'][()]
        sig=f['/filtered_{0:s}'.format(var)][()]
        sig_rms=f['/filtered_rms_{0:s}'.format(var)][()]
        print(len(tim))
        print(sig)
        f.close()

        for idx in enumerate(timearray):
            if idx%100==0:
                print('Reconstructing the signal at time {0:03.0f} of {1:03.0f}'.format(idx,len(timearray)))
            animated_base = Base()
            animated_base['0'] = Zone()
            animated_base[0].shared.connectivity = plane[0][0].connectivity
            orientation = 'spanwise' if np.std(plane[0][0]["x"]) < 1e-2 else 'streamwise'
            animated_base[0].shared["y"] = plane[0][0]["y"]
            animated_base[0].shared["z"] = plane[0][0]["z"]
            animated_base[0].shared["x"] = plane[0][0]["x"] 
            animated_base[0][str(idx)] = Instant()
            animated_base[0][str(idx)]['Filtered_'+var] = sig[:,idx]
            animated_base[0][str(idx)].attrs['Time'] = float(idx)
            w = Writer('hdf_antares')
            w['filename'] = os.path.join(target_dir,'Filtered_{0:s}_{1:03.0f}_{2:03.0f}_{3:03d}'.format(var,freq_min,freq_max,idx))
            w['base'] = animated_base
            w['dtype'] = 'float32'
            w.dump()
    text = 'Constructed FFT file is stored in {0:s}'.format(target_dir)
    print(f'\n{text}\n')  
    del w, plane
    text = 'Reconstruction Process Complete'
    print(f'\n{text:.^80}\n')  