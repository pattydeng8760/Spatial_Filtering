# Spatial_Filtering

A Python module for applying FFT-based filtering to transient CFD data (e.g. dilatation fields) extracted from LES or AVBP simulations. Designed for post-processing visualization.
Developed as part of advanced CFD/Aeroacoustics analysis workflows.

## Author

Patrick Deng
PhD Candidate – University of Toronto
Computational Fluid Dynamics & Aeroacoustics


## Features

- Extracts surface cutplanes from AVBP meshes using Antares
- Loads time-resolved HDF5 field data (e.g., `dilatation_node`)
- Applies FFT → Bandpass Filter → IFFT to isolate acoustic signatures
- Reconstructs filtered fields for visualization/analysis
- Modular design using Antares and NumPy stack


## Structure

spatial_filter/ 
    ├── init.py # Initializes module 
    ├── field_filter.py # Main entrypoint (FieldFilter class + CLI) 
    ├── mesh_utils.py # Surface extraction, cut mapping 
    ├── fft_utils.py # FFT and IFFT bandpass logic 
    ├── reconstruction.py # Reconstruction and export 
    ├── data_extraction.py # Temporal data loader

## How to Run

```bash
python spatial_filter/field_filter.py \
  --dir_to_post /path/to/postprocessed_h5_files \
  --mesh_dir /path/to/mesh_directory \
  --mesh_filename Bombardier_10AOA.mesh.h5 \
  --cut_location Cut_5mm_tip \
  --dt 3.1e-4 \
  --freq_min 800 \
  --freq_max 1500
```

## Outputs:
Filtered signals saved in ./Filtered_Signals/Cut_5mm_tip/
FFT logs saved to log_Cut_5mm_tip.txt

```bash

from spatial_filter.field_filter import FieldFilter
from spatial_filter import fft, create_tempo_base

args = argparse.Namespace(
    dir_to_post="...",
    mesh_dir="...",
    mesh_filename="...",
    cut_location="...",
    dt=3.1e-4,
    freq_min=800,
    freq_max=1500,
    zones=1500,
    override=False
)

filter_job = FieldFilter(args)
filter_job.run(args)
```

## Dependencies 
Antares (v1.18)
numpy
h5py
scipy
glob, os, sys, argparse

## Notes
Designed for structured/unstructured meshes from AVBP
Assumes Antares mesh extraction is compatible with cut locations
Frequencies are filtered using direct FFT masking for speed
Output HDF5 is Antares-readable (for ParaView visualization)

