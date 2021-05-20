# PhysiCell_intracellular: demonstrate Intracellular modeling, using SBML and libRoadrunner, in PhysiCell.

## Overview: 
PhysiCell is a flexible open source framework for building agent-based multicellular models in 3-D tissue environments.

**Reference:** A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellular Systems, PLoS Comput. Biol. 14(2): e1005991, 2018. DOI: [10.1371/journal.pcbi.1005991](https://dx.doi.org/10.1371/journal.pcbi.1005991)

Visit http://MathCancer.org/blog for the latest tutorials and help. 

## Intracellular demo model with libRoadrunner

Before you begin, you need to install the (binary) libRoadrunner. In the `/beta` directory, 
run `python setup_libroadrunner.py` and follow its instructions. We recommend
installing the (free) Anaconda Python distribution from https://www.anaconda.com/products/individual.

Next, try to run the following commands from a command line window.

```
$ make roadrunner_simple1
$ make
$ test_rr1
$ cd output
$ cp ../scripts/anim_svg_substrate_grid_step.py .
$ python anim_svg_substrate_grid_step.py 0 0 0 600 0 600 -300 300 -300 300 0
```
Rf. the `update_intracellular()` function in `main.cpp` to see what's happening between the PhysiCell substrate (oxygen) and the intracellular species (Oxy). This function is currently invoked from `main()` every `diffusion_dt`. Basically, that function just gets a cell's voxel value for oxygen, puts it into the SBML model, runs the Roadrunner solver, gets the new Oxy value, and copies it back into the voxel.
