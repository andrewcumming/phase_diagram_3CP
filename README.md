##Compiling

* The free-energy routines use cython to speed up the calculation. To compile  `free_energy.pyx` (2 component free energy) and `free_energy3.pyx`into modules, use

`python setup.py build_ext --inplace`
`python setup3.py build_ext --inplace`

To make an html to inspect the code you can use, e.g.
`cython -a free_energy.pyx`
which will make   `free_energy.html`

* `phase_diagram3_plots.py` uses [python-ternary](https://github.com/marcharper/python-ternary) to make ternary plots

* The code expects to find two directories `dat` and `out3`, so you should make those

`mkdir dat`

`mkdir out3`

##Running the code

`phase_diagram.py` calculates the phase diagram for a 2 component mixture. It generates a plot `phase_diagram_%d_%d.pdf` and the data is output in the `dat` directory

`phase_diagram3.py` calculates the phase diagram for a 3 component mixture. The data us output to the `dat` directory. 

`phase_diagram3_plots.py` produces ternary plots based on the data output from `phase_diagram3.py`. The plots are written to the `out3` directory
