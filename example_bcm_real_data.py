from __future__ import division
from bcm import *


cdata = cData(datafile='ingredlist',
              featurefile='ingrednames',
              namefile='recipe_names',
              is_binary_feature=True)

S = 5
F = len(cdata.x[0])
FG = [2] * F
N = len(cdata.x)
ITER = 200
s = iBCM(S, F, FG, N, ITER,
         prefix='recipe_example',
         alpha=1.0,
         feature_dictionary=cdata.feature_dic,
         data_name_dictionary=cdata.data_dic,
         q=0.4,
         is_saving_pickle=True,
         PEAKSIZE=20,
         lamb=2,
         is_debugging=False,
         is_verbose=True,
         is_visualize=True,
         mode='realdata_example')

if len(cdata.true_labels) > 0:
    s.sample_all(['z'], cdata)
else:
    s.sample_all([], cdata)
