Data and code for optimization of bond valence parameters in metalorganics
==========================================================================

Introduction
------------

This repository contains code and data needed to reproduce or tweak the large scale optimization of bond valence parameters for metalorganics taken from the Cambridge Structural Database.

Results based on this code have been published in: Zheng et al., *Data mining of Fe(II) and Fe(III) bond valence parameters, and their relevance for macromolecular crystallography*, Acta Crystallographica Section D, 2016.

There are two versions of the code:
* optimize-single.py: assumes a single parameter for each bond type
* optimize-split.py: allows fitting two parameters for bond types by splitting the input dynamically

Similarly, input and output data files are grouped into two directories: `data-single` and `data-split`.

Certain snapshots of selected files are imported into a [corresponding figshare article](http://figshare.com/articles/Data_and_code_for_optimization_of_bond_valence_parameters_in_metalorganics/964285), and each import is marked as a release/tag in this repository.
The DOI for the figshare article: [10.6084/m9.figshare.964285](http://dx.doi.org/10.6084/m9.figshare.964285)

How to rerun calculations
-------------------------

This repository relies on [a library that provides bond valence parameters](https://github.com/MinorLabUVa/bvparm)
which is linked as a submodule. We use these literature values as initial values for optimization.
To fetch this library, activate the submodule by issue the following command from this directory:
```bash
git submodule update
```

To run the single parameter optimization for a particular cation and oxidation state, run the following command:
```bash
python optimize-single.py <cation_name> <oxidation_state>
```

For example, to optimize the bond valence parameter for Na(I) based on the input binding sites:
```bash
python optimize-single.py sodium 1
```

Run the following command to optimize with two sets of parameters for Fe(II):
```bash
python optimize-split.py iron 2
```
We provide input data in this repository for optimizing with split parameters only for iron. The binding sites used for optimization can be filtered by the types of anions involved. In order to do this, adjust the variable `anion_selection` at the top of the `optimize-split.py` script.

Description of single parameter data files
------------------------------------------

In this case we provide input and output for all cations studied, which includes Ca(II), Fe(II), Fe(III), K(I), Na(I) and Zn(II). Both input and output files are in CSV format. The input file columns are:

1. Binding site ID
1. CSD accession code
1. The cation name
1. Assumed oxidation state of the cation
1. Coordination number
1. The name of the anion
1. Distance from cation to anion

For each binding site / accession code the number of rows will be the same as the coordinatin number.

The columns in the output files are:

1. CSD accession code
1. The cation name
1. The cation valence
1. The name of the anion
1. Distance from cation to anion
1. Valence attributed to this bond based on the optimize bond valence parameter.

Description of split parameter data
-----------------------------------

In this case we provide input and output only for Fe(II) and Fe(III), since these were the cations of interest and the only ones with significant changes when split parameters were used. Input files `bvparams2014d_iron2.csv` and `bvparams2014d_iron3.csv` contain the following columns:

1. Binding site accession code
1. The cation name
1. Distance from cation to anion
1. The anion name
1. Atomic numver of the anion

The output files, which match the pattern `optimized_iron<oxidation-state>_<anion_used>.csv`, have the following columns:

1. CSD accession code
1. Which bond valence parameter was chosen for this binding site (can be 0 or 1)
1. The cation name
1. The cation valence
1. Cumulative atomic valence (the total is the value for the last anion listed for the binding site)
1. The name of the anion
1. Distance from anion to cation
1. Optimized bond valence parameter for this bond
1. Valence attributed to this bond based on the optimized bond valence parameter

In addition to the above listings, the data files in this case include plots and their source data. Files following the pattern `plot_inital_iron<oxidation_state>_<anions used>` show the districution of bond valence sums based around the expected oxidation state based on literature parameters. Files following the patterm `plot_optimized_iron<oxidation_state>_<anions_used>` show the same distribution based on optimized parameters, and illustrate the convergence of the optimization proecdure. An example is shown below.

![Converge of two parameters for Fe-N binding sites and bond valence sum distribution](/data-split/plot_optimized_iron2_N_split.png "Optimization of two bond valence parameters for Fe-N homoleptic sites")
