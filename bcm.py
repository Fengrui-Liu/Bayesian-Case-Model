import sys
import os
from numpy import *
import copy as cp
import time
import itertools
from math import gamma, lgamma
import csv
import cPickle as pickle
import datetime
import argparse
from ctypes import cdll, c_int, c_double
import numpy.ctypeslib as npct
import json
import ctypes
import numpy.ctypeslib as npct
import glob
import operator
import platform
import matplotlib
import commands
from cData import *
import util
try:
    os.environ['DISPLAY']
    is_running_remotely = False
except:
    is_running_remotely = True
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

if 'Linux' in platform.system():
    shellfile_to_call = './compile_bcm_sharelib_ubuntu.sh'
    library_to_load = './libbcm.o'
elif 'Darwin' in platform.system():
    shellfile_to_call = './compile_bcm_sharelib.sh'
    library_to_load = './libbcm.dylib'
else:
    print 'Sorry, your OS is not supported.'
    sys.exit(1)

output = commands.getstatusoutput(shellfile_to_call)
if output[0] != 0:
    print 'Error in compiling bcm library: %s' % (output[1])
    sys.exit(1)
lib = cdll.LoadLibrary(library_to_load)

'''
  initialize all the functions in the lib
'''

# seed the randomness!
SEED = random.randint(0, sys.maxsize // 10000000000)
# for testing.
SEED = 1
random.seed(SEED)
lib.set_random_seed(c_int(SEED))

array_1d_double = npct.ndpointer(dtype=c_double, ndim=1, flags='CONTIGUOUS')
array_1d_int = npct.ndpointer(dtype=c_int, ndim=1, flags='CONTIGUOUS')
array_2d_int = npct.ndpointer(dtype=c_int, ndim=2, flags='CONTIGUOUS')
array_3d_int = npct.ndpointer(dtype=c_int, ndim=3, flags='CONTIGUOUS')
array_1d_double = npct.ndpointer(dtype=c_double, ndim=1, flags='CONTIGUOUS')
array_2d_double = npct.ndpointer(dtype=c_double, ndim=2, flags='CONTIGUOUS')
array_3d_double = npct.ndpointer(dtype=c_double, ndim=3, flags='CONTIGUOUS')

lib.normalize.restype = None
lib.normalize.argtypes = [array_1d_double, c_int, array_1d_double]
lib.get_al_val.restype = None
lib.get_al_val.argtypes = [
    c_int,
    c_int,
    c_double,
    c_double,
    c_int,
    array_1d_double]
lib.update_counting_variables.restype = None

lib.update_counting_variables.argtypes = [
    array_2d_int,
    array_3d_int,
    c_int,
    c_int,
    c_int,
    array_2d_int,
    array_2d_int,
    array_1d_int]

lib.sampleW_subroutine.restype = None
lib.sampleW_subroutine.argtypes = [
    c_int,
    c_int,
    c_int,
    array_2d_int,
    c_double,
    array_3d_int,
    array_2d_int,
    array_2d_int,
    c_double,
    c_int,
    array_1d_int,
    array_2d_int,
    array_2d_int,
    array_2d_double,
    array_3d_double,
    c_int]

lib.sampleZ_subroutine.restype = None
lib.sampleZ_subroutine.argtypes = [
    c_int,
    c_int,
    c_int,
    array_2d_int,
    c_double,
    array_3d_int,
    array_2d_int,
    array_2d_int,
    c_double,
    c_int,
    array_1d_int,
    array_2d_int,
    array_2d_int,
    array_3d_double,
    c_int]

lib.sampleP_subroutine.restype = None
lib.sampleP_subroutine.argtypes = [
    c_int,
    c_int,
    c_int,
    array_2d_int,
    c_double,
    array_3d_int,
    array_2d_int,
    array_1d_int,
    array_2d_int,
    c_double,
    c_int,
    array_1d_int,
    array_2d_int,
    array_2d_int,
    array_1d_double,
    array_2d_double,
    c_int]

lib.sample_all.restype = c_double
lib.sample_all.argtypes = [
    c_int,
    c_int,
    c_int,
    array_2d_int,
    c_double,
    array_3d_int,
    array_2d_int,
    array_1d_int,
    array_2d_int,
    c_double,
    c_int,
    array_1d_int,
    array_2d_int,
    array_2d_int,
    array_1d_double,
    array_2d_double,
    array_3d_double,
    array_2d_double,
    array_3d_double,
    c_int]

lib.sample_all_iterations.restype = c_double
lib.sample_all_iterations.argtypes = [
    c_int,
    c_int,
    c_int,
    array_2d_int,
    c_double,
    array_3d_int,
    array_2d_int,
    array_1d_int,
    array_2d_int,
    c_double,
    c_int,
    array_1d_int,
    array_2d_int,
    array_2d_int,
    array_1d_double,
    array_2d_double,
    c_int,
    array_3d_double,
    array_2d_double,
    array_3d_double,
    c_int]

lib.calculate_likelihood.restype = c_double
lib.calculate_likelihood.argtypes = [
    c_int,
    c_int,
    c_int,
    array_3d_int,
    array_2d_int,
    array_2d_int,
    c_double,
    c_int,
    array_1d_int,
    array_2d_double,
    array_2d_int,
    array_2d_int]

lib.get_scores.restype = None
lib.get_scores.argtypes = [
    c_int,
    c_int,
    c_int,
    array_2d_int,
    c_double,
    array_3d_int,
    array_2d_int,
    array_1d_int,
    array_2d_int,
    c_double,
    c_int,
    array_1d_int,
    array_2d_int,
    array_2d_int,
    array_1d_double,
    array_2d_double,
    array_3d_double,
    array_2d_double,
    array_3d_double]

lib.set_random_seed.restype = None
lib.set_random_seed.argtypes = [c_int]


class iBCM():

    def __init__(self, S, F, FG, N, ITERATION=300, alpha=0.1,
                 feature_dictionary=[],
                 feature_value_dictionary=[],
                 data_name_dictionary=[],
                 prefix='',
                 q=0.2,
                 PEAKSIZE=100,
                 lamb=2,
                 is_saving_pickle=True,
                 is_debugging=False,
                 is_verbose=True,
                 is_visualize=True,
                 mode='',
                 domain_name=''):
        self.S = S
        self.F = F
        if len(FG) != self.F:
            print 'number of feature values for each feature should be given as a list'
            sys.exit(1)

        # this should be a list of numbers for each feature
        self.FG = array(FG, dtype=int32)
        self.N = N

        self.p_val = zeros((self.S, self.F), dtype=int32)
        self.p_true_val = zeros((self.S, self.F), dtype=int32)
        self.p_index = zeros(self.S, dtype=int32)
        self.p_true_index = zeros(self.S, dtype=int32)
        self.w = ones((self.S, self.F), dtype=int32)
        self.w_true = ones((self.S, self.F), dtype=int32)

        self.phi = [[zeros(self.FG[f]) for f in xrange(self.F)]
                    for S in xrange(self.S)]
        self.phi_true = [[zeros(self.FG[f])
                          for f in xrange(self.F)] for S in xrange(self.S)]
        self.pi = zeros((self.N, self.S))
        self.pi_true = zeros((self.N, self.S))

        self.z = array([[random.randint(self.S) for f in xrange(self.F)]
                        for s in xrange(self.N)], dtype=int32)
        self.z_true = [] #array([[random.randint(self.S) for f in xrange(self.F)] for s in xrange(self.N)], dtype=int32)
        self.z_true_per_data = [] #array([random.randint(self.S) for s in xrange(self.N)], dtype=int32)

        self.x = zeros((self.N, self.F), dtype=int32)

        # Hyper params
        self.q = tile(q, (self.S, self.F))
        # good way to set this params for binary featured data might to set it to averge number of non zero elements of a
        # data devided by the entire number of data
        self.G0 = [1.] * self.N
        self.lamb_scalar = int32(lamb)
        self.lamb = [[self.lamb_scalar for fg in xrange(self.FG[f])] for f in xrange(
            self.F)]  # to select right alph, do self.lamb[f]
        self.alpha = float(alpha) / self.S
        # Peakier this is, easier for samplers to distinguish data.
        # defines how to peaky we want the phi to follow the prototype
        self.PEAKSIZE = int32(PEAKSIZE)
        self.likelihood = 0
        self.likelihood_array = []
        self.likelihood_handle = -1

        if len(feature_dictionary) > 0:
            if len(feature_dictionary.keys()) != self.F:
                print 'feature_dictionary keys (%s) must be the same length with the number of features (%s)' % (len(feature_dictionary.keys()), self.F)

                sys.exit(1)

            if len(feature_value_dictionary) > 0 and [ len( feature_value_dictionary[ feature_dictionary[f] ] ) for f in xrange(self.F) ] \
                    != self.FG:
                print 'feature_value_dictionary keys must be the same length with the number of feature values for each feature'
            self.feature_dictionary_inv = {
                feature_dictionary[key]: key for key in feature_dictionary.keys()}

        self.feature_dictionary = feature_dictionary
        self.feature_value_dictionary = feature_value_dictionary
        self.data_name_dictionary = data_name_dictionary

        ''' For Gibbs sampleing Stats
        '''
        self.nd = zeros(
            (self.N,
             self.S),
            dtype=int32)  # nd[i][k] = histogram(z[i]) within data count
        # x_by_sf[s][f][i] is not empty if and only if z[i][f] = s. If it is,
        # x_by_sf[s][f][i] = x[i][f] (int. fg val)
        self.nfg = array([[zeros(self.FG[f], dtype=int32) for s in xrange(
            self.S)] for f in xrange(self.F)])  # nfg[f][k] = hist( concat(x_by_sf[s][f]) )

        self.curr_iteration = 0
        self.ITERATION = ITERATION
        self.DEBUGGING = is_debugging
        self.is_verbose = is_verbose
        self.is_visualize = is_visualize
        self.data_for_debugging = []
        self.proto_set = [[] for s in xrange(self.S)]

        '''
          varialbes for explanations
        '''
        self.z_scores = array([[zeros(self.S, dtype=double) for f in xrange(
            self.F)] for n in xrange(self.N)], dtype=c_double)
        self.p_scores = zeros((self.S, self.N), dtype=c_double)
        self.w_scores = array([[zeros(2, dtype=double) for f in xrange(
            self.F)] for s in xrange(self.S)], dtype=c_double)

        if self.DEBUGGING:
            if self.is_verbose:
                print 'Generating simulation data...'
            self._generate_sim_data()

        self.is_saving_pickle = is_saving_pickle
        self.mode = mode
        '''
            visualize progress.
        '''
        self.fig = plt.figure(1, figsize=(20, 10))
        mode = 'sim' if self.DEBUGGING else 'real'
        ''' add some more info to file names '''
        self.postfix_for_files = datetime.date.today().strftime('%B%d') +\
            '_' + str(prefix) +\
            '_al-' + str(alpha) + '_iter-' + str(self.ITERATION) + '_nsubs-' + str(self.S) + \
            '_ndata-' + str(self.N) + '_nfeat-' + str(self.F) + '_nfv-' + str(self.FG[0]) +\
            '_q-' + str(self.q[0][0]) + '_PEAK-' + \
            str(self.PEAKSIZE) + '_lamb-' + str(self.lamb_scalar)

        self.results_folder = self.postfix_for_files
        os.system('mkdir ' + self.results_folder)
        protofile_name = self.results_folder + '/' + \
            'proto_' + self.postfix_for_files + '.txt'
        accfile_name = self.results_folder + '/' + \
            'acc_' + self.postfix_for_files + '.txt'
        if self.is_saving_pickle:
            picklefile_name = self.results_folder + '/' + \
                'pickle_' + self.postfix_for_files + '.txt'
            self.picklefile = open(picklefile_name, 'wb')

        self.accfile = open(accfile_name, 'w')
        self.protofile = open(protofile_name, 'w')
        self.protofile.write('SEED: ' + str(SEED) + '\n')

        self._make_all_nparray()

    def _w_true_exist(self):
        return True if len(self.w_true) > 0 else False

    def _z_true_exist(self):
        return True if len(self.z_true) > 0 else False

    def _z_true_per_data_exist(self):
        return True if len(self.z_true_per_data) > 0 else False

    def _p_true_val_exist(self):
        return True if len(self.p_true_val) > 0 else False

    def _p_true_exist(self):
        return True if len(self.p_true_index) > 0 else False

    def _phi_true_exist(self):
        return True if len(self.phi_true) > 0 else False

    def _generate_sim_data(self):
        '''
            genearte function generates a simulated data from BCM.
        '''

        if self.DEBUGGING:
            # or any([self.FG[f] < self.S for f in xrange(self.F)] ):
            if self.S > self.N:
                # TODO: is this still needed?
                print 'When debugging mode, the number of data should be bigger than the number of subspace'
                sys.exit(1)
            if True:
                self.data_for_debugging = util.get_data_for_debugging(
                    self.S,
                    self.FG,
                    self.F)
                if self.is_verbose:
                    print 'Printing data_for_debugging :%s' % (self.data_for_debugging)

        for s in xrange(self.S):
            self.w_true[s] = [
                random.binomial(
                    1,
                    self.q[s][f]) for f in xrange(
                    self.F)]
            if self.DEBUGGING:
                ''' This is for debugging the sampler. It generates diagonal pattern for w, so that
                    I can easily see if things are working well or not.
                '''
                #self.w_true[s] = [ 1 if (f%(self.S)==s or f%(self.S)==s+1 or  f%(self.S)==s+2 ) else 0  for f in xrange(self.F)]
                self.w_true[s] = [
                    1 if (
                        f == s or f == s +
                        1) else 0 for f in xrange(
                        self.F)]
                #self.w_true[s] = [ 1  for f in xrange(self.F)]

            ''' NOTE: tricky but possible situation to have all self.w[s] to be zero.
                This is case as if we don't have the subspace, so we can just remove them.
                Equivalently, here we make sure w has at least one zeros.
            '''
            while sum(self.w_true[s]) == 0:
                self.w_true[s] = [
                    random.binomial(
                        1,
                        self.q) for f in xrange(
                        self.F)]

            # At inference step, the data should be given, but for the generation step,
            # we aso generate this data point.
            self.p_true_val[s] = [
                random.randint(
                    0,
                    self.FG[f]) for f in xrange(
                    self.F)]
            if self.DEBUGGING:
                ''' This is for debugging the sampler. It generates diagonal pattern for w, so that
                    I can easily see if things are working well or not.
                    Assign the generated prototype to be my true prototype.
                '''
                self.p_true_val[s] = self.data_for_debugging[s]
                self.p_true_index[s] = s
                if self.is_verbose:
                    print 'p_true_val %s' % (self.p_true_val[s])
                    print 'w_true %s' % (self.w_true[s])
            for f in xrange(self.F):
                ''' normalize phi '''
                alval = zeros(self.FG[f], dtype=c_double)
                lib.get_al_val(
                    self.p_true_val[s][f],
                    self.w_true[s][f],
                    self.lamb_scalar,
                    self.PEAKSIZE,
                    self.FG[f],
                    alval)
                self.phi_true[s][f] = random.dirichlet(array(alval))

        self.z_true = array([[random.randint(self.S) for f in xrange(self.F)] for s in xrange(self.N)], dtype=int32)
        self.z_true_per_data =array([random.randint(self.S) for s in xrange(self.N)], dtype=int32)
        for n in xrange(self.N):
            self.pi_true[n] = random.dirichlet([self.alpha] * self.S)
            for f in xrange(self.F):
                self.z_true[n][f] = int32(
                    random.multinomial(
                        1,
                        self.pi_true[n]).argmax())
                self.x[n][f] = int32(
                    random.multinomial(
                        1, self.phi_true[
                            self.z_true[n][f]][f]).argmax())
            self.z_true_per_data[n] = max(
                set(self.z_true[n]), key=list(self.z_true[n]).count)

        if self.DEBUGGING:
            '''
                Now, since we generated the data by choosing data_for_debugging, we need to make sure those
                data points are actually included in the self.x.
                So we just insert them to some z_true and x[n]s.
            '''
            for s in xrange(self.S):
                self.z_true[s] = [s] * self.F
                self.x[s] = self.data_for_debugging[s]

                '''
                Find data that are generated that has the same important features as the prototypes
                They will basically generate the same score as the prototypes!
                '''
                imprtF = self.w_true[s].nonzero()
                self.proto_set[s] = [i for i in xrange(self.N) if all([a == b for a, b in zip(array(self.x[i])[imprtF],
                                                                                              array(self.data_for_debugging[s])[imprtF])]
                                                                      )]
        self._make_all_nparray()

    def _make_all_nparray(self):
        '''
        Make all into np array so that ctypes are happy when we hand it over to the cpp library!
        '''
        self.nd = array(self.nd, dtype=int32)
        self.nfg = array(self.nfg, dtype=int32)
        self.x = array(self.x, dtype=int32)
        self.z = array(self.z, dtype=int32)
        self.z_true = array(self.z_true, dtype=int32)
        self.p_true_val = array(self.p_true_val, dtype=int32)
        self.p_val = array(self.p_val, dtype=int32)
        self.p_true_index = array(self.p_true_index, dtype=int32)
        self.p_index = array(self.p_index, dtype=int32)
        self.w = array(self.w, dtype=int32)
        self.w_true = array(self.w_true, dtype=int32)
        self.G0 = array(self.G0, dtype=c_double)
        self.q = array(self.q, dtype=c_double)
        self.z_scores = array(self.z_scores, dtype=c_double)
        self.p_scores = array(self.p_scores, dtype=c_double)
        self.w_scores = array(self.w_scores, dtype=c_double)

    def _initializeVariables(self, data=[], addNoiseVar=['w', 'p', 'z']):
        '''
         initializeVariables funtion is called AFTER latent variables are randomly initialized.
         This is to check if there is ground truth, add noise to the ground truth if requested
         and update the counting variables for inference (nfg and nd).
        '''
        if not isinstance(addNoiseVar, list):
            print 'the third input to initializeVariables should be list of strings of variables'
            sys.exit(1)

        if self.DEBUGGING == True:
            self.z = array(cp.deepcopy(self.z_true), dtype=int32)
            self.p_val = cp.deepcopy(self.p_true_val)
            self.p_index = cp.deepcopy(self.p_true_index)
        self.w = cp.deepcopy(self.w_true)

        if isinstance(data, cData) and len(
                data.x) > 0:  # working with real data
            if len(data.z_true) == 0 and 'z' in addNoiseVar:
                print 'z_true does not exist! cannot add noise to z'
                sys.exit(1)
            if len(data.w_true) == 0 and 'w' in addNoiseVar:
                print 'w_true does not exist! cannot add noise to w'
                sys.exit(1)
            if len(data.p_true_val) == 0 and 'p' in addNoiseVar:
                print 'p_true does not exist! cannot add noise to p'
                sys.exit(1)
            if len(data.z_true) > 0:
                self.z_true = array(data.z_true, dtype=int32)
                self.z = array(self.z_true, dtype=int32)
            if len(data.z_true_per_data) > 0:
                self.z_true_per_data = array(data.z_true_per_data, dtype=int32)
            if len(data.w_true) > 0:
                self.w_true = array(data.w_true)
                self.w = self.w_true
            if len(data.p_true_val) > 0:
                self.p_true_val = array(data.p_true_val)
                self.p_val = self.p_true_val
            if len(data.p_true_index) > 0:
                self.p_true_index = array(data.p_true_index)
                self.p_index = self.p_true_index
            self.x = array(data.x, dtype=c_int)
        '''
          if DOHALFMIX is true, the gibbs samples starts from more favorable location.
        '''
        DOHALFMIX = False
        if 'w' in addNoiseVar:
            if DOHALFMIX:
                BADSTART_IND = [0, 1, self.S - 1]
                for i in BADSTART_IND:
                    self.w[i] = [1 for f in xrange(self.F)]
            else:
                self.w = [[1 for f in xrange(self.F)] for s in xrange(self.S)]
        if 'p' in addNoiseVar:
            if DOHALFMIX:
                BADSTART_IND = [0, 1, self.S - 1]
                for i in BADSTART_IND:
                    idchoosen = random.randint(self.N)
                    self.p_index[i] = idchoosen
                    self.p_val[i] = self.x[idchoosen]
            else:
                self.p_index = [random.randint(self.N) for s in xrange(self.S)]
                self.p_val = [self.x[self.p_index[s]] for s in xrange(self.S)]
        if 'z' in addNoiseVar:
            if DOHALFMIX:
                # every FREQ_BADSTART data points we will randomly shuffle z.
                FREQ_BADSTART = 2
                changeRows = array(
                    [1 if i % FREQ_BADSTART == 0 else 0 for i in xrange(self.N)]).nonzero()[0]
                if self.is_verbose:
                    print 'Adding noise to z. Every %s is noise ' % (FREQ_BADSTART)
                self.z[changeRows] = [
                    [
                        random.randint(
                            self.S) for f in xrange(
                            self.F)] for i in xrange(
                        len(changeRows))]
                ''' then, for the rows that are not completely changed, add noise too. '''
                FREQ_BADSTART_2 = 10
                self.z = array([[random.randint(self.S) if f % FREQ_BADSTART_2 == 0 else self.z[
                               i][f] for f in xrange(self.F)] for i in xrange(len(self.z))], dtype=int32)
            else:
                self.z = array(
                    [[random.randint(self.S) for f in xrange(self.F)] for i in xrange(self.N)])
        if 'z' not in addNoiseVar:
            self.phi = cp.deepcopy(self.phi_true)
            #self.z = array([ [ random.randint(self.S) for f in xrange(self.F) ] for i in xrange(self.N) ])

        self._make_all_nparray()
        self._update_counting_variables()

    def _update_counting_variables(self):
        lib.update_counting_variables(
            self.nd,
            self.nfg,
            self.N,
            self.S,
            self.F,
            self.z,
            self.x,
            self.FG)

    @staticmethod
    def compute_phi_standalone(lamb, PEAKSIZE, p, w, S, F, FG, nfg):
        phi = [[zeros(FG[f]) for f in xrange(F)] for s in xrange(S)]
        for s in xrange(S):
            for f in xrange(F):
                al = zeros(FG[f], dtype=c_double)
                lib.get_al_val(p[s][f], w[s][f], lamb, PEAKSIZE, FG[f], al)
                phi[s][f] = (array(nfg[f][s]) + al) / \
                    (float(sum(nfg[f][s])) + sum(array(al)))
                if (sum(phi[s][f]) - 1) > 1e-10:
                    print 'compute phi error: sum of phi[s][f] needs to be 1'
                    sys.exit()
        return array(phi)

    def getNd(self, i, s):
        return self.nd[i][s]

    def getNfg(self, f, s, fg):
        return self.nfg[f][s][fg]

    def sampleZ(self):
        '''
           NOTE: only sampleing Z depends on PEAKSIZE and how distinguishable each prototype is
                making data_for_debugging 000, 111 makes eaiser for z to sample.
                also, if PEAKSIZE is big, it also make sampling Z easy.
                having small number of S helps
                having all features turned on (Self.q[s][f] = 1) helps.
        '''
        lib.sampleZ_subroutine(self.N, self.S, self.F, self.nd, self.alpha,
                               self.nfg, self.p_val, self.w, self.lamb_scalar, self.PEAKSIZE, self.FG,
                               self.x, self.z, self.z_scores, SEED)

    def sampleP(self):
        '''
          sample prototytpes
        '''
        lib.sampleP_subroutine(self.N, self.S, self.F, self.nd, self.alpha,
                               self.nfg, self.p_val, self.p_index, self.w, self.lamb_scalar,
                               self.PEAKSIZE, self.FG,
                               self.x, self.z, self.G0, self.p_scores, SEED)

    def sampleW(self):
        '''
          sample subspaces
        '''
        lib.sampleW_subroutine(self.N, self.S, self.F, self.nd, self.alpha, self.nfg,
                               self.p_val, self.w, self.lamb_scalar, self.PEAKSIZE, self.FG, self.x, self.z, self.q, self.w_scores, SEED)

    def sample_all(self, addnoiseVars, data=[]):

        if not isinstance(data, cData):
            data = self.x
        self._calculate_likelihood()
        if self.is_verbose:
            print 'initializing variables'
        self._initializeVariables(data, addnoiseVars)
        if self.is_verbose:
            print 'DONE initializing variables'

        if not is_running_remotely and self.is_visualize:
            self.fig.show()
            if self.is_verbose:
                print 'showing figure'

        if self.is_saving_pickle:
            pickle.dump(self.x, self.picklefile)
        self._calculate_likelihood()
        SAMPLEALL=False
        if not SAMPLEALL:
            for i in xrange(self.ITERATION):
                self.curr_iteration = i
                if (i % 5 == 0 and i > 0):
                    accs = self.print_save_acc(i)
                    if self.is_visualize:
                      self.visualize(accs, i)
                self.write_prototypes(i)
                if self.is_verbose:
                    print '. iteration %s .' % (i)
                self.sampleP()
                self.sampleW()
                self.sampleZ()
                self._calculate_likelihood()

                lib.get_scores(self.N, self.S, self.F, self.nd, self.alpha,
                               self.nfg, self.p_val, self.p_index, self.w,
                               self.lamb_scalar, self.PEAKSIZE, self.FG,
                               self.x, self.z, self.G0, self.q, self.z_scores, self.p_scores, self.w_scores)

        else:
            # if we do not need visualization, do all the iterations from cpp
            # side. It's faster.
            self._make_all_nparray()
            if self.is_verbose:
                print 'sampling all iterations from the cpp library... (nothing will print until done)'
            lib.sample_all_iterations(self.N, self.S, self.F, self.nd, self.alpha,
                                      self.nfg, self.p_val, self.p_index, self.w, self.lamb_scalar, self.PEAKSIZE, self.FG,
                                      self.x, self.z, self.G0, self.q, self.ITERATION, self.z_scores, self.p_scores, self.w_scores, SEED)
            if self.is_verbose:
                print 'done sampling %s times' % (self.ITERATION)
            self.print_save_acc(-1)
        self.accfile.close()
        self.protofile.close()

        if self.is_saving_pickle:
            self._calculate_likelihood()
            self.save_picklefile_and_close()

    def save_picklefile_and_close(self):
        pickle.dump({'p_index': self.p_index, 'p_val': self.p_val,
                     'z': self.z, 'w': self.w, 'nd': self.nd, 'nfg': self.nfg,
                     'lamb_scalar': self.lamb_scalar, 'PEAKSIZE': self.PEAKSIZE,
                     'N': self.N, 'S': self.S, 'F': self.F,
                     'FG': self.FG, 'q': self.q,
                     'pi': self.pi, 'phi': self.phi,
                     'p_scores': self.p_scores, 'w_scores': self.w_scores, 'z_scores': self.z_scores,
                     'likelihood': self.likelihood,
                     'p_true_val': self.p_true_val,
                     'p_true_index': self.p_true_index,
                     'w_true': self.w_true,
                     'phi_true': self.phi_true,
                     'pi_true': self.pi_true,
                     'z_true': self.z_true,
                     'z_true_per_data': self.z_true_per_data,
                     'data_for_debugging': self.data_for_debugging,
                     }, self.picklefile)
        self.picklefile.close()

    def _calculate_likelihood(self):
        '''
            Mainly to show that Gibbs is moving to the right direction: try on sim data first!
        '''
        self.likelihood = lib.calculate_likelihood(
            self.N,
            self.S,
            self.F,
            self.nfg,
            self.p_val,
            self.w,
            self.lamb_scalar,
            self.PEAKSIZE,
            self.FG,
            self.q,
            self.z,
            self.x)
        self.likelihood_array.append(self.likelihood)

    @staticmethod
    def compute_pi_standalone(N, S, z):
        pi = zeros((N, S))
        for n in xrange(N):
            pi_t = array([list(z[n]).count(i) for i in xrange(S)])
            pi[n] = pi_t / float(sum(pi_t))
        return pi

    def print_save_acc(self, iter):
        '''
          print (if verbose mode) and save accuracy and likelihood
        '''
        if self.DEBUGGING:

            self.accfile.write('iter: ' + str(iter))
            wsacc = ''
            phiacc = ''
            pacc = ''
            zacc = ''

            wsacc = sum(sum(array([[aa == bb for aa, bb in zip(
                self.w[s], self.w_true[s])] for s in xrange(self.S)]))) / float(self.S * self.F)
            self.accfile.write(' ws ' + str('%.3f' % (wsacc)))
            phiacc = sum(array([[self.phi[s][f].argmax() == self.phi_true[s][f].argmax(
            ) for s in xrange(self.S)] for f in xrange(self.F)]).flatten()) / float(self.S * self.F)
            self.accfile.write(' phi ' + str('%.3f' % (phiacc)))

            p_est_remap = [
                s if self.p_index[s] in self.proto_set[s] else self.p_index[s] for s in xrange(
                    self.S)]
            pacc = sum([pp == pp_true for pp,
                        pp_true in zip(p_est_remap,
                                       self.p_true_index)]) / float(self.S)
            self.accfile.write(' p ' + str('%.3f' % (pacc)))
            zacc = (sum(sum(array([[a == b for a, b in zip(self.z[i], self.z_true[
                    i])] for i in xrange(self.N)]))) / float(self.N) / float(self.F))
            self.accfile.write('  z ' + str('%.3f' % (zacc)))
            self.accfile.write(' likelihood:' +
                               str('%.4f' %
                                   (self.likelihood)) +
                               '\n')

            if self.is_verbose:
                print '--Accuracy-- p: %.3f w: %.3f phi: %.3f z: %.2f likelihood: %.2f' % (pacc, wsacc, phiacc, zacc, self.likelihood)
            return {'w': wsacc, 'phi': phiacc, 'p': pacc, 'z': zacc}

        elif self._z_true_exist():
            self.accfile.write('iter: ' + str(iter))
            zacc = (sum(sum(array([[a == b for a, b in zip(self.z[i], self.z_true[
                    i])] for i in xrange(self.N)]))) / float(self.N) / float(self.F))
            self.accfile.write('  z ' + str('%.3f' % (zacc)))
            self.accfile.write(' likelihood:' +
                               str('%.4f' %
                                   (self.likelihood)) +
                               '\n')
            if self.is_verbose:
                print '--Accuracy-- z: %.2f likelihood: %.2f' % (zacc, self.likelihood)
            return {'w': '', 'phi': '', 'p': '', 'z': zacc}
        else:
            if len(set(self.z_true_per_data)) > 2:
                acc = util.get_acc(
                    self.z,
                    self.p_index,
                    self.N,
                    self.S,
                    self.z_true_per_data,
                    self.F)
                if self.is_verbose:
                    print '--Accuracy-- acc:%.2f likelihood: %.2f' % (acc, self.likelihood)
                self.accfile.write('unsuper ' +
                                   str('%.3f' %
                                       (acc)) +
                                   ' likelihood ' +
                                   str(self.likelihood) +
                                   '\n')
            else:
                self.accfile.write('iter: ' + str(iter))
                self.accfile.write(' likelihood ' +
                                   str('%.4f' %
                                       (self.likelihood)) +
                                   '\n')
            # just returning meaningless thing.
            return {'w': '', 'phi': '', 'p': '', 'z': ''}

    def write_prototypes(self, i):
        '''
          print and write the current prototypes into a file
        '''
        self.protofile.write(
            ' -------------------------------- iter ' +
            str(i) +
            ' ---------------------------------\n')
        ''' if metafeature_dictionary exists, it is using meta features! '''
        if len(self.feature_dictionary) > 0:
            for s in xrange(self.S):
                if len(self.data_name_dictionary) > 0:
                    proto_in_words = [
                        self.data_name_dictionary[
                            self.p_index[s]]]  # if self.p_index[s] else ' ' ]
                else:
                    proto_in_words = str(self.p_index[s]) + ' data point'
                data_feature_indices = array(self.p_val[s]).nonzero()[0]
                feature_vals_protoType = [
                    self.feature_dictionary[i] for i in data_feature_indices]
                subspace_in_words = [
                    self.feature_dictionary[i] for i in array(
                        self.w[s]).nonzero()[0] if i in data_feature_indices] if sum(
                    self.w[s]) != 0 else ' None '
                if len(data_feature_indices) > 0:
                    ratio_subspace = float(
                        len(subspace_in_words)) / float(len(data_feature_indices))
                else:
                    ratio_subspace = 'NA'

                if self._z_true_exist():
                    add_label_str = ' label ' + \
                        str(self.z_true_per_data[self.p_index[s]])
                else:
                    add_label_str = ''

                str_to_write_protofile = 'subspace ' + str(s) + '  prototypes :======================= ' +  \
                    str(proto_in_words) + add_label_str + ' subspace ratio: ' + \
                    str(ratio_subspace) + ' ============================\n'

                self.protofile.write(str_to_write_protofile)
                self.protofile.write(' subspace in words ' + str(subspace_in_words) + '\n' +
                                     ' Feature Vals of prototype \n' + '\n'.join(feature_vals_protoType) +
                                     '\n')

    def visualize(self, accs, curr_iter):
        '''
            visualize current samples.
        '''
        if self.DEBUGGING:
            subplot_row = 4
            subplot_col = 2

            plt.clf()
            title_for_figure = self.postfix_for_files[:len(
                self.postfix_for_files) / 2] + '\n' + self.postfix_for_files[len(self.postfix_for_files) / 2:]
            self.fig.suptitle(
                'Currently storing the results at: ' +
                title_for_figure,
                fontsize=10)

            if self._p_true_exist():
                max_width_of_p_pretty = max(self.p_true_index) + 1
                p_true_binary = [[1 if i == self.p_true_index[s] else 0 for i in xrange(
                    max_width_of_p_pretty)] for s in xrange(self.S)]
                #p_est_binary = [ [ 1 if i==self.p_index[s] else 0 for i in xrange(self.N) ] for s in xrange(self.S) ]
                a_p_true = self.fig.add_subplot(subplot_row, subplot_col, 1)
                a_p_true.set_title('true p')
                a_p_true.imshow(
                    p_true_binary,
                    interpolation='nearest',
                    cmap=cm.Greys_r)
                a_p_true.axis('off')

            ''' Here map all p to be protoset '''
            p_est_remap = [
                s if self.p_index[s] in self.proto_set[s] else self.p_index[s] for s in xrange(
                    self.S)]
            max_width_of_p_pretty = self.S  # self.N
            p_est_binary = [[1 if i == p_est_remap[s] else 0 for i in xrange(
                max_width_of_p_pretty)] for s in xrange(self.S)]
            a_p_est = self.fig.add_subplot(subplot_row, subplot_col, 2)
            a_p_est.set_title('est. p')
            a_p_est.imshow(
                p_est_binary,
                interpolation='nearest',
                cmap=cm.Greys_r)
            a_p_est.axis('off')

            if self._w_true_exist():
                a_w_true = self.fig.add_subplot(subplot_row, subplot_col, 3)
                a_w_true.set_title('true w')
                a_w_true.axis('off')
                implot_w_true = a_w_true.imshow(
                    array(
                        self.w_true),
                    interpolation='nearest',
                    cmap=cm.Greys_r)
                implot_w_true.set_clim(0., 1.)
            a_w_est = self.fig.add_subplot(subplot_row, subplot_col, 4)
            a_w_est.set_title('est. w')
            a_w_est.axis('off')
            implot_w_est = a_w_est.imshow(
                array(
                    self.w),
                interpolation='nearest',
                cmap=cm.Greys_r)
            implot_w_est.set_clim(0., 1.)
            # print self.w_true
            # print self.w
            if self._z_true_exist():
                a_z_true = self.fig.add_subplot(subplot_row, subplot_col, 5)
                a_z_true.set_title('true z')
                a_z_true.axis('off')
                a_z_true.imshow(array(self.z_true).reshape(-1,
                                                           max(self.F,
                                                               self.N)),
                                interpolation='nearest',
                                cmap=cm.Greys_r)

            a_z_est = self.fig.add_subplot(subplot_row, subplot_col, 6)
            a_z_est.set_title('est. z')
            a_z_est.imshow(array(self.z).reshape(-1,
                                                 max(self.F,
                                                     self.N)),
                           interpolation='nearest',
                           cmap=cm.Greys_r)
            a_z_est.axis('off')

            if len(self.data_for_debugging) > 0:
                a_data_for_debugging_true = self.fig.add_subplot(
                    subplot_row,
                    subplot_col,
                    7)
                a_data_for_debugging_true.set_title('data for debugging')
                a_data_for_debugging_true.imshow(
                    self.data_for_debugging,
                    interpolation='nearest',
                    cmap=cm.Greys_r)
                a_data_for_debugging_true.axis('off')

            a_likelihood = self.fig.add_subplot(subplot_row, subplot_col, 8)
            a_likelihood.set_title('log likelihood')

            a_likelihood.plot(
                range(len(self.likelihood_array)), self.likelihood_array)
            # add accuracy string in the figure.
            self.fig.text(0, 0, 'w:' + str('%.3f ' % (accs['w'])) + ' p:' + str('%.3f' % (accs['p'])) +
                          ' z:' + str('%.3f' % (accs['z'])) + ' phi:' + str('%.3f' % (accs['phi'])) + ' likelihood:' + str('%.3f' % (self.likelihood)))
            plt.draw()
        else: #if len(self.z_true) > 0:
            Dict = self.get_dict_of_currn_vars()
            fig = self.fig
            '''
              visualize p and w
          '''
            subplot_row = 4 #+ self.S // 2 + 1
            if self._z_true_exist():
              subplot_col = 2
            else:
              subplot_col = 1

            plt.clf()
            title_for_figure = self.postfix_for_files[:len(
                self.postfix_for_files) / 2] + '\n' + self.postfix_for_files[len(self.postfix_for_files) / 2:]
            fig.suptitle(title_for_figure, fontsize=10)
            first_fig_index = 1
            a_p_est = fig.add_subplot(subplot_row, subplot_col, first_fig_index)
            a_p_est.set_title('est. p')
            a_p_est.imshow([[i if self.p_index[s] == i else 0 for i in xrange(
                self.N)] for s in xrange(self.S)], interpolation='nearest', cmap=cm.Greys_r)
            a_p_est.axis('off')

            a_w_est = fig.add_subplot(subplot_row, subplot_col, first_fig_index + subplot_col)
            a_w_est.set_title('est. w')
            a_w_est.axis('off')
            implot_w_est = a_w_est.imshow(
                array(
                    self.w),
                interpolation='nearest',
                cmap=cm.Greys_r)
            implot_w_est.set_clim(0., 1.)
            a_z_est = fig.add_subplot(subplot_row, subplot_col, first_fig_index + 2*subplot_col)
            a_z_est.set_title('est. z')
            a_z_est.imshow(array(self.z).reshape(-1,
                                                 max(self.F,
                                                     self.N)),
                           interpolation='nearest',
                           cmap=cm.Greys_r)
            a_z_est.axis('off')

            a_likelihood = fig.add_subplot(subplot_row, subplot_col, first_fig_index + 3*subplot_col)
            a_likelihood.set_title('log likelihood')

            a_likelihood.plot(
                range(len(self.likelihood_array)), self.likelihood_array)
            # add accuracy
            fig.text(0, 0, 'w:' + str(accs['w']) + ' p:' + str(accs['p']) +
                     ' z:' + str(accs['z']) + ' phi:' + str(accs['phi']) + ' likelihood:' + str(self.likelihood))
            if self._z_true_exist():
              a_z_true = fig.add_subplot(subplot_row, subplot_col, first_fig_index + 2*subplot_col + 1)
              a_z_true.set_title('true z')
              a_z_true.axis('off')
              a_z_true.imshow(array(self.z_true).reshape(-1,
                                                        max(self.F,
                                                            self.N)),
                              interpolation='nearest',
                              cmap=cm.Greys_r)
            plt.draw()
        '''
        else:
            plt.clf()
            subplot_row = 1
            subplot_col = 1

            title_for_figure = self.postfix_for_files[:len(
                self.postfix_for_files) / 2] + '\n' + self.postfix_for_files[len(self.postfix_for_files) / 2:]
            self.fig.suptitle(title_for_figure, fontsize=10)
            max_width_of_p_pretty = max(self.p_true_index) + 1
            a_likelihood = self.fig.add_subplot(subplot_row, subplot_col, 0)
            a_likelihood.set_title('log likelihood')
            a_likelihood.plot(
                range(len(self.likelihood_array)), self.likelihood_array)
            plt.draw()
            self.fig.text(0, 0, 'iter ' +
                          str(curr_iter) +
                          ' likelihood:' +
                          str(self.likelihood))
        '''

        if self.fig:
            self.fig.savefig(
                self.results_folder +
                '/Gibbs_' +
                str(curr_iter) +
                'iter.png',
                dpi=80)

    def get_dict_of_currn_vars(self):
        '''
        returns a dictionary of the current prototype and subspaces.
        out = {'prototype_name':{ 'subspaces': [0,1, 0], # subspaces
                    'cluster_id': cluster_id,
                    'features_in_strings': ['green', 'circle','roundpattern'],  # features in string
                    'feature_indices': [feature id] #
                    }
              }
        '''
        proto_dict = {}
        alldata = set()
        # if true, it shows the max pscored prototype, not the one that is
        # sampled.
        # update pi
        self.pi = iBCM.compute_pi_standalone(self.N, self.S, self.z)
        self.phi = array(
            self.compute_phi_standalone(
                self.lamb_scalar,
                self.PEAKSIZE,
                self.p_val,
                self.w,
                self.S,
                self.F,
                self.FG,
                self.nfg),
            dtype=c_double)
        SHOW_P_INDEX = True
        if SHOW_P_INDEX:
            p_indices = self.p_index
        else:
            p_indices = [self.p_scores[cluster_id].argmax()
                         for cluster_id in xrange(self.S)]

        for cluster_id, proto_index in enumerate(p_indices):
            if len(self.data_name_dictionary) > 0:
                proto_name = self.data_name_dictionary[
                    proto_index].replace("'", "").lower()
            else:
                proto_name = str(proto_index)

            subspaces = list(
                self.w[cluster_id][
                    self.x[proto_index].nonzero()[0]])
            feature_strings = [
                self.feature_dictionary[k].replace(
                    "'",
                    "") for k in self.x[proto_index].nonzero()[0]]
            feature_ids = [self.feature_dictionary_inv[fs]
                           for fs in feature_strings]

            proto_dict[cluster_id] = {
                'prototype': proto_name,
                'subspaces': subspaces,
                'features_in_strings': feature_strings,
                'feature_indices': feature_ids}

        return proto_dict
