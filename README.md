## Simulation code associated with the article "Complete Singularity Analysis for the Perspective-Four-Point Problem" submitted to IJCV

### Dependencies:
- [`ViSP`](https://visp.inria.fr) (available as Ubuntu package libvisp-dev)
- Shipped in when calling `git clone https://github.com/oKermorgant/singularity_4points --recurse-submodules`: 
  - [`log2plot`](https://github.com/oKermorgant/log2plot) to create figures
  - [`lambdatwist_pnp`](https://github.com/oKermorgant/lambdatwist_pnp.git) for the P4P library
  - [`opengv`](https://github.com/laurentkneip/opengv) for UPnP and EPnP
  
### Executables
- rose: the rose pattern simulation
- sphere: the growing sphere simulation
- visual_servo: the multi-start visual servoing simulation
- energy: the energy landscape computation

All executables are configured from the `config.yaml` file. Singularity scenes are from 1 to 10. Non-singular scenes are built from scenes 1 and 2 and are named 1o and 2o.

If the `dataPath` folder does not exists, result files will be created in `<source directory>/results`

### Scripts
In the `analyze` folder, some scripts allow calling multiple simulations at once. 

In particular, `spheres.sh` will run the `sphere` simulation for all pose computation methods + all refinement methods. This may take quite some time depending on the computer.

The Python scripts are to be run after the simulations are over, in order to create summary plots.
