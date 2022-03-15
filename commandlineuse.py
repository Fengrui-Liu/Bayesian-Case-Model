from __future__ import division

from bcm import *

if __name__=='__main__':
    '''
    To run animal data, do:
          python bcm.py -file animals.csv animals_features.csv animals_names.csv animals_labels.csv
    '''
    ''' data '''
    parser = argparse.ArgumentParser(description="processing BCM options")
    parser.add_argument('-mode', default='simple', required=True, choices=['sim','real'], help='simple/datafile')
    parser.add_argument('-S', type=int, default=4, help='number of clusters')
    parser.add_argument('-F',  type=int, default=4, help='number of features')
    parser.add_argument('-FG',  type=int, default=4, help='number of feature values per feature')
    parser.add_argument('-N',  type=int, default=100, help='number of data points')

    ''' model params '''
    parser.add_argument('-ITER',  type=int, default=100, help='number of iterations')
    parser.add_argument('-alpha',  type=float, default=0.1, help='Dirichlet parameter')
    parser.add_argument('-PEAKSIZE',  type=float, default=10, help='Dirichlet parameter')
    parser.add_argument('-lamb',  type=float, default=1, help='Dirichlet parameter')
    parser.add_argument('-q', type=float, default=0.8, help='Bernoulli parameter')

    ''' data files '''
    parser.add_argument('-file', default='',  help='data file')
    parser.add_argument('-featurefile', default='',  help='file contains names of features')
    parser.add_argument('-namefile', default='',  help='file contains names of data points')
    parser.add_argument('-labelfile', default='',  help='file contains labels of data points')
    parser.add_argument('-resultFile_prefix', default=None,  help='file contains labels of data points')
    parser.add_argument('-savepickle', default=True,  help='True only if you want to save pickle files (big)')
    parser.add_argument('-picklefileToLoad', default='',  help='a pickle file to load data from')
    parser.add_argument('-subsample', default=0,  help='whether to subsample the dataset to have less features')
    ''' other options '''
    parser.add_argument('-visualize', default=1, choices=['1','0'],  help='visualize current Gibbs samples')
    parser.add_argument('-verbose', default=1,choices=['1','0'],  help='print useful stuff')

    args = parser.parse_args()

    if args.mode == 'datafile' and args.file == '':
      print 'provide the file -file FILE' 
      sys.exit(1)

    args.resultFile_prefix = args.mode if not args.resultFile_prefix else args.resultFile_prefix

    if 'sim' == args.mode:
        #tic = time.clock()
        s = iBCM(args.S, args.F, [args.FG]*args.F, args.N, args.ITER, \
            prefix = args.resultFile_prefix, \
            alpha = args.alpha, q = args.q, 
            is_saving_pickle =args.savepickle,\
            PEAKSIZE = args.PEAKSIZE,\
            lamb = args.lamb,\
            is_debugging = True,\
            is_verbose=bool(int(args.verbose)),\
            is_visualize=bool(int(args.visualize)),\
            mode=args.mode)
        DOINFERENCE = True 
        if DOINFERENCE:
          s.sample_all(['p','w','z'])
          #s.sample_all(['z'])
    elif 'real' == args.mode: #datafile provided
      doSubsample=False
      if 'animals' in args.file and args.subsample == 1:
        doSubsample=True

      cdata = cData(args.file, \
            featurefile = args.featurefile,\
            namefile = args.namefile,\
            labelfile = args.labelfile,\
            is_binary_feature = True)
      if args.labelfile:
        S = len(set(cdata.true_labels))
      else:
        S = args.S
 
      F = len(cdata.x[0])
      FG = [2]*F
      N = len(cdata.x)
      s = iBCM(S, F, FG, N, args.ITER, \
          prefix = args.resultFile_prefix, \
          alpha = args.alpha, \
          feature_dictionary = cdata.feature_dic, \
          data_name_dictionary=cdata.data_dic, \
          q = args.q, \
          is_saving_pickle = args.savepickle, \
          PEAKSIZE = args.PEAKSIZE,\
          lamb = args.lamb,\
          is_debugging = False,\
          is_verbose = bool(int(args.verbose)),\
          is_visualize=bool(int(args.visualize)),\
          mode=args.mode)

      if len(cdata.true_labels) > 0:
        s.sample_all(['z'], cdata)
      else:
        s.sample_all([], cdata)


