from __future__ import division
import itertools
import numpy as np

def get_data_for_debugging(S, FG, F):
    # fg  < s: then each s becomes one fg.
    if FG[0] > S:
      return np.array([[ s for f in xrange(F) ] for s in xrange(S) ], dtype=np.int32)
    else:
      # find permutations of m number of elements for a subset of fv, such
      # that the number of possible permunations are greater than s.
      permutations = []
      m = 0
      while len(permutations) < S:
        len_fg = len(range(FG[0]))
        if m > len_fg:
          array_to_permute = range(FG[0]) + list(range(FG[0])*(m-len_fg))[:(m-len_fg)] 
        else:
          array_to_permute = range(FG[0])
        permutations = list(set(list(itertools.permutations(array_to_permute, m))))
        m += 1
      # repleat permutations[s] F//len_permute times
      return np.array([ list(permutations[s]*(F//len(permutations[s])+1))[:F]  for s in xrange(S) ], dtype=np.int32)


def get_acc(z, p_index, numX, numS, true_labels, F):
    subspace_labels =  true_labels[ p_index ]
    # just duplicate the true label..
    true_labels_allFeatures = np.array( [  [t]*F for t in true_labels] ).flatten()
    estimated_labels = np.array( [ subspace_labels[ [ zif for zif in z[i]] ] for i in xrange(numX) ] ).flatten()
    print 'subspace true labels %s' % (subspace_labels)
    print 'subspace p_index %s ' % (p_index)
    acc = sum( estimated_labels == true_labels_allFeatures)/float(len(estimated_labels))
    return acc




