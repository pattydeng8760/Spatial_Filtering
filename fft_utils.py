# fft_utils.py
import numpy as np
import os
import h5py
import scipy.fftpack as sf
import scipy.signal as sg
from multiprocessing import Pool, cpu_count
import multiprocessing as mp


# functions to find the right factors to segmet the FFT data
# Helper functions used in the code
def find_factor(number):
    """
    Find and return a list of factors for a given number.
    """
    return [i for i in range(1, number + 1) if number % i == 0]

def closest(lst, K):
    """
    Return the element in lst that is closest to K.
    """
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]
# The FFT function
def fft(output, nnodes, var_interest, dt, freq_min, freq_max, zones, override=False):
    """ 
    Function to perform FFT on the selected data and save the filtered data between the selected frequency range
    Args:
        output (str): The output directory where the FFT data is stored
        nnodes (int): The number of nodes in the data
        var_interest (list): The list of variables to be processed
        dt (float): The time step size
        freq_min (float): The minimum frequency for filtering
        freq_max (float): The maximum frequency for filtering
        zones (int): The number of zones for FFT processing
        override (bool, False): If True, overrides the number of probes. Recommended to be False.
    Returns:
        outfile (str): The path to the output file containing the filtered data
    """
    text = 'Beginning the FFT process'
    print(f'\n{text:.^60}\n')  
    factor = find_factor(nnodes)
    if override == True:
        nprobe = int(np.ceil(nnodes/zones))
    else:
        nprobe = int(np.ceil(nnodes/zones))
        nprobe = int(closest(factor,nprobe))
    fmin=np.array((freq_min,))
    fmax=np.array((freq_max,))
    print(f'Filtering the selected data from {freq_min :3d} Hz to {freq_max :3d}  Hz \n\n') 
    if not os.path.isdir(output):
        sys.error('Output directory does not exist.')
    # Loop over the variables of interest
    for var in var_interest:
        print('Compute FFT for variable {0:s}'.format(var))
        filename='{0:s}/field_{1:s}.h5'.format(output,var)  
        f=h5py.File(filename,'r')
        data=f['/field'][()]
        f.close()
        nb_ite=data.shape[1]
        nnodes=data.shape[0]
        nlen=nb_ite if nb_ite%2==0 else nb_ite+1
        fft_len=nlen
        npart=int(nnodes/nprobe)
        full_fft_vec=np.zeros((nnodes,fft_len),dtype='complex64')
        freq=sf.fftfreq(nlen,d=dt)[0:fft_len]
        scale_factor=1.0/float(nlen)
        print(f'    The number of total nodes in the plane is: {nnodes :.1f}.') 
        print(f'    The number of timesteps investigated is: {nlen :.0f}.') 
        print(f'    The number of points in a region is: {nprobe :.1f}.') 
        print(f'    The number of discrete sub-domains in the plane is: {npart :.1f}.\n') 
        if npart*nprobe == nnodes:
            print('    The total number of nodes is equal to the number of subdomain * number of nodes in each subdomain:')
            print(f'        {npart :3d} * {nprobe :3d} = {nnodes :3d}\n\n')
        else:
            print('\n\n Warning!!')
            print('\n\n The number the number of subdomain * number of nodes in each subdomain is NOT equal to the total number of nodes')
            print(f'{npart :3d} * {nprobe :3d} = {nprobe*nprobe :3d} =/= {nnodes :3d}\n\n')
        #============================   
        # Forward FFT procedure
        #============================   
        print("Beginning the forward FFT proecdure.")
        for ipart in range(npart):
            if ipart == npart-1 or ipart % 100 == 0:
                print('    Computing fft for mesh part {0:d}'.format(ipart))
            istart=ipart*nprobe
            iend=(ipart+1)*nprobe
            p_fluct=sg.detrend(data[istart:iend,:],axis=1,type='constant')
            fft_vec=sf.fft(p_fluct,axis=1)
            full_fft_vec[istart:iend,:]=fft_vec[:,0:fft_len]
        istart=npart*nprobe
        iend=nnodes
        p_fluct=sg.detrend(data[istart:iend,:],axis=1,type='constant')     
        fft_vec=sf.fft(p_fluct,axis=1)
        full_fft_vec[istart:iend,:]=fft_vec[:,0:fft_len]
        full_fft_vec=full_fft_vec*scale_factor
        print('Forward FFT Complete')
        
        #============================   
        # Inverse FFT procedure with bandpass filter
        #============================  
        print("\nStarting the inverse-FFT procedure.")
        # IFFT procedure: loop over the frequency ranges
        for cur_freq_min, cur_freq_max in zip(fmin, fmax):
            # Allocate space for the arrays
            filtered_fft_vec = np.zeros((nnodes, fft_len), dtype='complex64')
            filtered_ifft_vec = np.zeros((nnodes, fft_len), dtype='complex64')
            filtered_ifft_vec_rms = np.zeros((nnodes), dtype='complex64')
            timesignal = np.linspace(0, dt * nlen, nlen)

            # Determine the frequency indices within the desired range
            freq_ind = ((freq >= cur_freq_min) * (freq <= cur_freq_max)) | ((freq >= -cur_freq_max) * (freq <= -cur_freq_min))
            filtered_fft_vec[:, freq_ind] = full_fft_vec[:, freq_ind]

            # Process each subdomain (mesh part)
            for ipart in range(npart):
                if ipart == npart-1 or ipart % 100 == 0:
                    print('    Computing ifft for mesh part {0:d}'.format(ipart))
                istart = ipart * nprobe
                iend = (ipart + 1) * nprobe

                temp_ifft_vec = sf.ifft(filtered_fft_vec[istart:iend, :])
                filtered_ifft_vec[istart:iend, :] = temp_ifft_vec[:, 0:fft_len]
                filtered_ifft_vec_rms[istart:iend] = np.std(temp_ifft_vec, axis=1)

            # After processing all subdomains, save the results to file
            outfile = '{0:s}/new_ifft_{1:s}_freq_{2:03.0f}_{3:03.0f}.h5'.format(output, var, cur_freq_min, cur_freq_max)
            with h5py.File(outfile, 'w') as g:
                g['filtered_{0:s}'.format(var)] = np.real(filtered_ifft_vec.real)
                g['filtered_rms_{0:s}'.format(var)] = np.real(filtered_ifft_vec_rms.real)
                g['time'] = timesignal
            print('Inverse-FFT Complete')
    text = 'Constructed FFT file is stored in {0:s}'.format(outfile)
    print(f'\n{text}\n')  
    text = 'FFT Process Complete'
    print(f'\n{text:.^60}\n')  
    return outfile

def segment_fft_worker(data_segment, dt, seg_id):
    """
    Worker function to compute FFT on a pre-sliced data segment.
    
    Args:
        data_segment (ndarray): Pre-sliced data array for this segment.
        dt (float): Time step size.
        seg_id (int): Segment ID (for logging).
    Returns:
        ndarray: The FFT result of the segment.
    """
    if seg_id % 100 == 0:
        print(f"    Computing FFT for segment {seg_id}")
    p_fluct = sg.detrend(data_segment, axis=1, type='constant')
    return np.real(sf.fft(p_fluct, axis=1)),seg_id

def segment_ifft_worker(fft_segment, seg_id):
    """
    Worker function to compute IFFT on a filtered FFT segment.
    
    Args:
        fft_segment (ndarray): Pre-sliced FFT data for this segment.
        seg_id (int): Segment ID (for logging).
    Returns:
        ndarray: The IFFT result of the segment.
    """
    if seg_id % 100 == 0:
        print(f"    Computing IFFT for segment {seg_id}")
    return np.real(sf.ifft(fft_segment, axis=1)),seg_id


def parallel_fft(output, nnodes, var_interest, dt, freq_min, freq_max, zones,  num_cpu=10):
    """
    Performs FFT-based filtering in parallel by partitioning the node data,
    applying FFT and IFFT on each segment, and then stitching the results.
    
    Args:
        output (str): Directory containing the HDF5 field data and where to save output.
        nnodes (int): Total number of nodes in the field.
        var_interest (list): List of variable names (e.g., ['dilatation_node']) to process.
        dt (float): Time step size.
        freq_min (float): Minimum frequency for the bandpass filter.
        freq_max (float): Maximum frequency for the bandpass filter.
        zones (int): Approximate number of segments (blocks) to partition the node data.
        num_cpu (int): Number of CPU cores to use.
    
    Returns:
        outfile (str): The path to the output file with the filtered results.
    """
    # Use the spawn context for multiprocessing
    text = 'Beginning the FFT filtering process'
    print(f'\n{text:.^60}\n') 
    num_cpu = mp.cpu_count() if num_cpu is None else num_cpu
    print(f"Using {num_cpu} CPU cores for parallel processing.")
    # For each variable, perform the parallel FFT/IFFT filtering.
    for var in var_interest:
        print(f"Processing variable: {var}")
        outfile = os.path.join(output, f'ifft_{var}_freq_{int(freq_min):03}_{int(freq_max):03}.h5')
        # Check if the output directory exists
        if os.path.exists(outfile):
            print('Output fft data already exists. Skipping processing...')
        else:
            print(f"\n---->Extracting field data for variable: {var}")
            filepath = os.path.join(output, f'field_{var}.h5')
            # Open the file once in the parent process to determine shape and pre-slice the data.
            with h5py.File(filepath, 'r') as f:
                dataset_raw = np.array(f['/field'])
                data_shape = dataset_raw.shape  # (nnodes, ntime)
                # Confirm that nnodes matches the provided value (optional)
                if data_shape[0] != nnodes:
                    print(f"Warning: nnodes ({nnodes}) does not match dataset shape ({data_shape[0]})")
                ntime = data_shape[1]
                # Ensure an even number of time steps for FFT
                fft_len = ntime if ntime % 2 == 0 else ntime - 1
            dataset = np.array(dataset_raw,dtype=np.float32)
            del dataset_raw  # Free up memory
            
            # Pre-slice the data along the node dimension.
            node_indices = np.array_split(np.arange(nnodes), zones)
            # Pre-read each slice from the file (only the slice, not the whole dataset).
            fft_tasks = [(dataset[indices, :fft_len],dt,block_num) for block_num,indices in enumerate(node_indices)]
            print(f"Pre-sliced {nnodes} nodes into {len(fft_tasks)} segments.")
            
            # Frequency and scaling setup.
            freq = sf.fftfreq(fft_len, d=dt)
            scale_factor = 1.0 / fft_len
            freq_ind = ((freq >= freq_min) & (freq <= freq_max)) | ((freq >= -freq_max) & (freq <= -freq_min))
            timesignal = np.linspace(0, dt * fft_len, fft_len)
            
            # Allocating memory for the dataset
            full_fft_vec = np.zeros((nnodes, len(freq)), dtype=np.float32)
            full_ifft_vec = np.zeros((nnodes, fft_len), dtype=np.float32)
            
            #=================================
            # Parallel forward FFT: process each pre-sliced segment.
            #=================================
            print("\n---->Computing forward FFT ...")
            with Pool(processes=num_cpu) as pool:
                fft_results = pool.starmap(segment_fft_worker, fft_tasks)
            # Stitch FFT segments together.
            for fft_results_block, seg_id in fft_results:
                # Store the FFT results in the full array.
                full_fft_vec[node_indices[seg_id], :] = fft_results_block
            full_fft_vec *= scale_factor  # Apply scaling factor to the full FFT vector
            # Apply bandpass filter: keep only the desired frequency range.
            filtered_fft_vec = np.zeros((nnodes, len(freq)), dtype=np.float32)
            filtered_fft_vec[:, freq_ind] = full_fft_vec[:, freq_ind]
            del full_fft_vec  # Free up memory
            
            # Prepare IFFT tasks: split the filtered FFT vector in the same segments.
            ifft_tasks = [(filtered_fft_vec[indices, :], seg_id) for seg_id, indices in enumerate(node_indices)]
            
            
            #=================================
            # Parallel inverse FFT: process each filtered segment.
            #=================================
            print("\n---->Computing inverse FFT ...")
            with Pool(processes=num_cpu) as pool:
                ifft_results = pool.starmap(segment_ifft_worker, ifft_tasks)
            # Stitch IFFT segments together.
            for ifft_results_block, seg_id in ifft_results:
                # Store the IFFT results in the full array.
                full_ifft_vec[node_indices[seg_id], :] = ifft_results_block
                
            # Save the final filtered (IFFT) result and RMS to an output HDF5 file.
            with h5py.File(outfile, 'w') as g:
                g[f'filtered_{var}'] = full_ifft_vec
                #g[f'filtered_rms_{var}'] = filtered_ifft_vec_rms.real
                g['time'] = timesignal
            print(f"\n---->Filtered data for variable '{var}' saved to:\n {outfile}")
    
    text = 'FFT filtering process Complete'
    print(f'\n{text:.^60}\n')  
    return outfile

def butter_filter_worker(data_segment, dt, fmin, fmax, order, seg_id):
    """
    Worker function to apply a Butterworth bandpass filter on a data segment.
    Uses filtfilt for zero-phase filtering.
    """
    b, a = sg.butter(order, [fmin / (0.5/dt), fmax / (0.5/dt)], btype="band")
    if seg_id % 100 == 0:
        print(f"    Filtering segment {seg_id}")
    filtered = sg.filtfilt(b, a, data_segment, axis=1)
    return filtered, seg_id
        
        
def parallel_butterworth(output, nnodes, var_interest, dt, freq_min, freq_max, zones, num_cpu=5):
    """
    Performs Butterworth bandpass filtering in parallel by partitioning the node data,
    applying the filter on each segment, and then stitching the results.
    
    Args:
        output (str): Directory containing the HDF5 field data and where to save output.
        nnodes (int): Total number of nodes in the field.
        var_interest (list): List of variable names (e.g., ['dilatation_node']) to process.
        dt (float): Time step size.
        freq_min (float): Lower cutoff frequency.
        freq_max (float): Upper cutoff frequency.
        zones (int): Approximate number of segments (blocks) to partition the node data.
        order (int): Order of the Butterworth filter.
        num_cpu (int): Number of CPU cores to use.
    
    Returns:
        outfile (str): The path to the output file with the filtered results.
    """
    text = "Beginning the Butterworth filtering process"
    print(f"\n{text.center(80, '.')} \n")
    order = 4
    for var in var_interest:
        print(f"Processing variable: {var}")
        filepath = os.path.join(output, f"field_{var}.h5")
        outfile = os.path.join(output, f"butter_{var}_order_{order}_fmin_{int(freq_min):03}_fmax_{int(freq_max):03}.h5")
        if os.path.exists(outfile):
            print("Output filtered data already exists. Skipping processing...")
            continue
        print(f"\n----> Extracting field data for variable: {var}")
        # Open file and load data into memory (adjust if memory is an issue)
        with h5py.File(filepath, "r") as f:
            dataset = np.array(f["/field"])
            data_shape = dataset.shape  # (nnodes, ntime)
            if data_shape[0] != nnodes:
                print(f"Warning: nnodes ({nnodes}) does not match dataset shape ({data_shape[0]})")
            ntime = data_shape[1]
        # Pre-slice the data along the node dimension.
        node_indices = np.array_split(np.arange(nnodes), zones)
        butter_tasks = [(dataset[indices, :], dt, freq_min, freq_max, order, block_num)
            for block_num, indices in enumerate(node_indices)]
        print(f"Pre-sliced {nnodes} nodes into {len(butter_tasks)} segments.")
        
        # Design Butterworth filter parameters.
        nyq = 0.5 / dt
        low = freq_min / nyq
        high = freq_max / nyq
        
        print("----> Parallelizing Butterworth filtering ...")
        with Pool(processes=num_cpu) as pool:
            results = pool.starmap(butter_filter_worker, butter_tasks)
        # Stitch filtered segments together.
        filtered_data = np.zeros((nnodes, ntime), dtype=np.float32)
        for seg_result, seg_id in results:
            filtered_data[node_indices[seg_id], :] = seg_result
        # Optionally, compute RMS (standard deviation over time) for each node.
        rms = np.std(filtered_data, axis=1)
        timesignal = np.linspace(0, dt * ntime, ntime)
        
        # Save the filtered results and RMS to an output HDF5 file.
        with h5py.File(outfile, "w") as g:
            g[f"filtered_{var}"] = filtered_data
            g[f"filtered_rms_{var}"] = rms
            g["time"] = timesignal
        print(f"----> Filtered data for variable '{var}' saved to:\n {outfile}")
    
    text = "Butterworth filtering process complete"
    print(f"\n{text.center(80, '.')} \n")
    return outfile
