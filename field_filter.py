# field_filter.py
import argparse
import os
import sys
import builtins
from mesh_utils import extract_surface
from data_extraction import create_tempo_base
from fft_utils import fft
from reconstruction import reconstruction

class FieldFilter:
    """A class to handle the filtering of field data from AVBP cutplanes.
    Args:
        dir_to_post (str): The directory where the AVBP cutplane files are located.
        mesh_dir (str): The directory where the AVBP mesh file is located.
        mesh_filename (str): The name of the AVBP mesh file.
        cut_location (str): The location of the cut plane, defined either explicitly via PIV planes
            or with a %_TE location designating the distance from the trailing edge.
        output (str): The destination directory where the extracted data is to be stored (mesh, baseline data, post-fft data).
        target_dir (str): The directory where the filtered data will be saved.
        var_interest (list): The list of variables to be extracted from the AVBP cutplane files.
        dt (float): The time step for the data extraction.
        freq_min (float): The minimum frequency for filtering.
        freq_max (float): The maximum frequency for filtering.
        zones (int): The number of zones for FFT processing.
        override (bool): If True, overrides the number of probes. Recommended to be False.
    """
    def __init__(self, args):
        self.dir_to_post = args.dir_to_post
        self.mesh_dir = args.mesh_dir
        self.mesh_filename = args.mesh_filename
        self.cut_location = args.cut_location
        self.output = f'{args.cut_location}'
        self.freq_min = args.freq_min
        self.freq_max = args.freq_max
        self.target_dir = os.path.join(self.output, args.cut_location+'_Filtered_'+str(freq_min)+'_'+str(freq_max))
        self.var_interest = ['dilatation_node']
        self.dt = args.dt
        self.zones = args.zones
        self.override = args.override
        os.makedirs(self.target_dir, exist_ok=True)
        # Redirecting stdout to a log file
        log_file = f'log_{args.cut_location}.txt'
        sys.stdout = open(log_file, "w", buffering=1)
        def print_redirect(text): builtins.print(text); os.fsync(sys.stdout)
        self.print = print_redirect

    def run(self):
        self.print(f'\n{"Starting the Filtering Program":.^120}\n')
        mesh, mesh_nodes = extract_surface(self.mesh_dir, self.mesh_filename, self.output, 'EXTRACT', self.cut_location)
        nnodes, nb_times = extract_data(self.dir_to_post, self.cut_location, self.output, self.var_interest, reload_data=False)
        fft_file = fft(self.cut_location, self.output, nnodes, self.var_interest, self.dt, self.freq_min, self.freq_max, self.zones, self.override)
        reconstruction(self.dir_to_post, mesh, self.var_interest, nb_times, self.freq_min, self.freq_max, fft_file, self.target_dir)
        self.print(f'\n{"Filtering Program Complete":.^120}\n')

def parse_args():
    parser = argparse.ArgumentParser(description='Run FFT filtering pipeline on transient data.')
    parser.add_argument('--dir_to_post', required=True, help='Directory containing .h5 time-resolved fields')
    parser.add_argument('--mesh_dir', required=True, help='Directory of AVBP mesh')
    parser.add_argument('--mesh_filename', required=True, help='Name of AVBP mesh file (e.g. Bombardier_10AOA.mesh.h5)')
    parser.add_argument('--cut_location', required=True, help='Name of cut location')
    parser.add_argument('--dt', type=float, required=True, help='Time step between snapshots')
    parser.add_argument('--var_interest', nargs='+', default=['dilatation_node'], help='Variables of interest to extract')
    parser.add_argument('--freq_min', type=int, default=800, help='Min frequency for bandpass')
    parser.add_argument('--freq_max', type=int, default=1500, help='Max frequency for bandpass')
    parser.add_argument('--zones', type=int, default=1500, help='Number of spatial segments in FFT')
    parser.add_argument('--override', type=bool, default=False, help='Force override of FFT segmentation zones')
    return parser.parse_args()

def main():
    args = parse_args()
    field_filter = FieldFilter(args)
    field_filter.run()

if __name__ == '__main__':
    main()
    sys.stdout.close()
    sys.stdout = sys.__stdout__
    print("Filtering process completed. Check the log file for details.")

