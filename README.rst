This is a backup of Bayesian Case Model Python Package from `Duke Interpretable ML Lab <https://users.cs.duke.edu/~cynthia/papers.html>`_ while the original repository was missing.


# Bayesian Case Model Python Package (original doc) #
======================================================

This is a python implementation of Bayesian Case Model, a clustering model that can learn prototypical examples and subspaces (i.e., important features) to explain each of the learned clusters. For details, please refer to: http://people.csail.mit.edu/beenkim/KimRudinShahNIPS2014

Currently only supporting Mac and Linux.

## Installing using easy_install ##

easy_install pybcm

## Manual build ##

Clone the repository:

```bash

git clone --recursive git@github.com:Fengrui-Liu/Bayesian-Case-Model.git
cd pybcm
python setup.py install --with-cython

```


## Running BCM ##

There are two modes -- simulation (-mode sim) and real data (-mode real).

# Command line #

For simulation,

```bash
python commandlineuse.py -mode sim -S 2 -F 2 -FG 2 -N 100 -ITER 400
```
where S is the number of clusters, F is the number of features, FG is the types a feature can take (e.g., 2 for binary features), N is the number of data to generate and ITER is number of Gibbs sampling iterations.

For real data,

```bash
python commandlineuse.py -mode real -file DATAFILE -featurefile FEATURENAME_FILE -namefile DATANAME_FILE -ITER 100 -S 3
```

where DATAFILE is the data file, FEATURENAME_FILE is a file where each row corresponds to the name of the feature (optional) and DATANAME_FILE is a file that contains names of the data points (optional). For the exact formats of these files, please refer to the example directory.

To view the full options, do

```bash
python commandlineuse.py -h
```
## Dataset ##

We also provide the recipe dataset used in the paper.

ingredlist: a comma separated values indicates the existence of each feature for each data point.
ingrednames: each row indicates the names of each feature
recipenames: each row indicates names of each recipe (data points)

To run the recipe dataset:

```bash
python commandlineuse.py -mode real -file ingredlist -featurefile ingrednames -namefile recipenames -ITER 1000 -S 5
```


## References ##

* Been Kim, Cynthia Rudin and Julie Shah [The Bayesian Case Model: A Generative Approach for Case-Based Reasoning and Prototype Classification]


@inproceedings{kim2014bayesian,
  title={The {B}ayesian {C}ase {M}odel: A Generative Approach for Case-Based Reasoning and Prototype Classification},
  author={Kim, Been and Rudin, Cynthia and Shah, Julie A},
  booktitle={Advances in Neural Information Processing Systems},
  pages={1952--1960},
  year={2014}
}

## Author ##

[Been Kim] (https://github.com/Beenie)
