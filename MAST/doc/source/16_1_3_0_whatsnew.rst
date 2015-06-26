#############################
What's new in version 1.3.0
#############################

**Additions:**

* STEM image processing additions to structopt. See **NO TAG YET**
* Charge-density-based pathfinding for NEB. See **NO TAG YET**
* Optional structure indexing (atomic indexing) for structures which change. Se :doc:`3_1_1_structure`.
* 14-frequency concentrated diffusion model in the Diffusion Coefficient post-processing tool. See ** ???**

**Fixes:**


**Changes for users:**

* Scaling is its own section in the input file, and is no longer under the structure section. See :doc:`3_1_9_scaling` for format changes.
* Instead of a mast_recipe.log file in each recipe directory, there is now one large mast.log file in $MAST_CONTROL.
* Instead of recipe-level trapped errors causing MAST to fail, MAST will create a MAST_ERROR document in the recipe directory and log a warning to ``$MAST_CONTROL/mast.log``

**Changes for programmers:**

* In order to see some logged messages that are now logged to the DEBUG level, set environment variable MAST_DEBUG to any value.
