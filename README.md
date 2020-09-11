Calculates the phase diagram for a three-component classical plasma. This code uses analytic fits to free energies of the liquid and solid phase to construct the free energy of the mixture. The double tangent construction is used to identify unstable regions in the phase diagram.

### Compiling

* The free-energy routines use cython to speed up the calculation. To compile  `free_energy.pyx` (2 component free energy) and `free_energy3.pyx` (3 component free energy) into modules, use

`python setup.py build_ext --inplace`

`python setup3.py build_ext --inplace`

To make an html to inspect the code you can use, e.g.
`cython -a free_energy.pyx`
which will make   `free_energy.html`

* `phase_diagram3_plots.py` uses [python-ternary](https://github.com/marcharper/python-ternary) to make ternary plots

* The code expects to find two directories `dat` and `out3`, so you should make those

`mkdir dat`

`mkdir out3`

### Running the code

`phase_diagram.py` calculates the phase diagram for a 2 component mixture. It generates a plot `phase_diagram_%d_%d.pdf` and the data is output in the `dat` directory

`phase_diagram3.py` calculates the phase diagram for a 3 component mixture. The data us output to the `dat` directory. 

`phase_diagram3_plots.py` produces ternary plots based on the data output from `phase_diagram3.py`. The plots are written to the `out3` directory

### References

* [Medin & Cumming 2010](https://ui.adsabs.harvard.edu/abs/2010PhRvE..81c6107M/abstract)

* [Caplan et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...860..148C/abstract)
