#!/bin/sh
set -e
echo "======================================================================"
echo " RUNNING CLEAN-ROOM ISOLATED VERIFICATION SUITE (python3 -I)"
echo "======================================================================"
echo "\n[1/4] Verifying Exhibit I: Renormalization-Free Lattice QFT..."
/usr/bin/python3 -I Run_1_standalone_ultrafinitist_qft.py --quick --out-dir Run_1_plots

echo "\n[2/4] Verifying Exhibit II: Ultrafinitist Classical Mechanics Over Q..."
/usr/bin/python3 -I Run_2_standalone_rational_mechanics.py --steps 1000 --plots-dir Run_2_plots

echo "\n[3/4] Verifying Exhibit III Symbolic Spatial & No-Blowup Bounds Over Q..."
/usr/bin/python3 -I Run_3_verify_spatial_and_bounds.py

echo "\n[4/4] Verifying Exhibit III: Singularity-Free 2D/3D Fluid Turbulence Over Q..."
/usr/bin/python3 -I Run_3_standalone_ultrafinitist_turbulence.py

echo "\n======================================================================"
echo " ALL THREE ULTRAFINITIST EXHIBITS VERIFIED CLEANLY (EXIT 0)"
echo "======================================================================"
