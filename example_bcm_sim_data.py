from __future__ import division
from bcm import *

S = 3  # number of clusters.
F = 3  # number of features
FG = 2  # number of values a feature can take (2 for binary features)
N = 300  # number of data to generate
ITER = 300  # number of Gibbs iterations
resultFile_prefix = 'simulation_example'
alpha = 0.1
q = 0.4
lamb = 1
PEAKSIZE = 10
is_verbose = True
is_visualize = True

s = iBCM(S, F, [FG] * F, N, ITER,
         prefix=resultFile_prefix,
         alpha=alpha, q=q,
         is_saving_pickle=True,
         PEAKSIZE=PEAKSIZE,
         lamb=lamb,
         is_debugging=True,
         is_verbose=is_verbose,
         is_visualize=is_visualize,
         mode='example')

s.sample_all(['p', 'w', 'z'])
