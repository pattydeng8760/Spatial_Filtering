# fft_utils.py
import numpy as np
import h5py
import scipy.fftpack as sf
import scipy.signal as sg

# functions to find the right factors to segmet the FFT data
def find_factor(number):
    """Find the factors of a number
    """
    factor = []
    count = 0
    for i in range(1,number+1):
        if number % i == 0:
            factor = np.append(factor,i)
            count += 1
    return factor
            
def closest(lst, K):
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]

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
    print(f'\n{text:.^80}\n')  
    factor = find_factor(nnodes)
    if override == True:
        nprobe = int(np.ceil(nnodes/zones))
    else:
        nprobe = int(np.ceil(nnodes/zones))
        nprobe = int(closest(factor,nprobe))
    fmin=np.array((freq_min,))
    fmax=np.array((freq_max,))
    print(f'Filtering the selected data from {freq_min :.2f} Hz to {freq_max :.2f}  Hz \n\n') 
    if not os.path.isdir(output):
        sys.error('Output directory does not exist.')
    # Loop over the variables of interest
    for var in var_interest:
        print('Compute FFT for variable {0:s} in {1:s}'.format(var,nmp))
        filename='{0:s}/field_{1:s}.h5'.format(output,var)  
        f=h5py.File(filename,'r')
        data=f['/field'][()]
        f.close()
        nb_ite=data.shape[1]
        nnodes=data.shape[0]
        nlen=nb_ite if nb_ite%2==0 else nb_ite+1
        
        print("Beginning the forward FFT proecdure.")
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
            print('\n\n    The domain segment is sucessful!')
            print('    The total number of nodes is equal to the number of subdomain * number of nodes in each subdomain.')
            print(f'        {npart :3d}*{nprobe :3d}={nnodes :3d}\n\n')
        else:
            print('\n\n Warning!!')
            print('\n\n The number the number of subdomain * number of nodes in each subdomain is NOT equal to the total number of nodes')
            print(f'{npart :3d}*{nprobe :3d}={nprobe*nprobe :3d}=/={nnodes :3d}\n\n')
        #============================   
        # Forward FFT procedure
        #============================   
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
        print("\n\nStarting the inverse-FFT procedure.")
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
    print(f'\n{text:.^80}\n')  
    return outfile