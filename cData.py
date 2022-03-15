from __future__ import division
import csv
from numpy import *


class cData():

    def __init__(self, datafile, origfile=None, featurefile=None,
                 metafeaturefile=None,
                 namefile=None,
                 labelfile=None,
                 make_dataname_dict=False,
                 is_binary_feature=True,
                 z_true=[],
                 w_true=[], p_true_val=[], p_true_index=[], z_true_per_data=[]):

        data = cData.convert_2_int(cData.readCSV(datafile))
        if origfile:
            self.origdata = cData.convert_2_int(cData.readCSV(origfile))
            self.origdata_feature_index = cData.get_index_none_zeros(
                self.origdata)
        else:
            self.origdata = []
            self.origdata_feature_index = []
        featurenames = cData.make1D(cData.remove_quotes_and_lower(
            cData.remove_trailing_spaces(cData.readCSV(featurefile))))
        self.feature_dic = cData.make_dic(featurenames)
        if metafeaturefile:
            self.meta_feature_to_orig_mapping = cData.get_metafeature_to_orig_feature_mapping(
                cData.readCSV(metafeaturefile))
        else:
            self.meta_feature_to_orig_mapping = []

        self.true_labels = []
        if labelfile:
            self.true_labels = cData.make1D(
                cData.convert_2_int(cData.readCSV(labelfile)))

        if len(data) == 0:
            print ' no data '
            sys.exit(1)
        if is_binary_feature and len(set(array(data).flatten())) > 2:
            self.x, self.x_feature_order = cData.get_feature_order_and_data(
                data)
        elif is_binary_feature and len(set(array(self.origdata).flatten())) > 2:
            self.x, self.x_feature_order = cData.get_feature_order_and_data(
                data, origdata=self.origdata)
        else:
            self.x = data
            self.x_feature_order = []

        self.data_dic = {}
        if namefile:
            self.data_dic = cData.make_dic(cData.make1D(cData.remove_quotes_and_lower(
                cData.remove_trailing_spaces(cData.readCSV(namefile)))))
        elif make_dataname_dict:
            self.data_dic = cData.make_dic(
                [str(i) for i in xrange(len(self.x))])

        F = len(self.x[0])
        self.z_true = [
            [tt] * F for tt in self.true_labels] if len(self.true_labels) > 0 else []
        self.z_true_per_data = self.true_labels if len(
            self.true_labels) > 0 else []
        self.w_true = w_true if len(w_true) > 0 else []
        self.p_true_val = p_true_val if len(p_true_val) > 0 else []
        self.p_true_index = p_true_index if len(p_true_index) > 0 else []
        #self.feature_dictionary = feature_dictionary if len(feature_dictionary) >0 else []
        #self.featureValueDictionary = featureValueDictionary if len(featureValueDictionary) >0 else []

    @staticmethod
    def get_feature_order_and_data(x, origdata=None):
        if origdata:
            x_feature_order = [[x_tt for x_tt in x_t if x_tt > 0]
                               for x_t in origdata]
        else:
            x_feature_order = [[x_tt for x_tt in x_t if x_tt > 0] for x_t in x]
        data = [[1 if x_tt > 0 else 0 for x_tt in x_t] for x_t in x]
        return data, x_feature_order

    @staticmethod
    def readCSV(filename):
        out = []
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) > 0:
                    out.append([rr for rr in row if len(rr) > 0])
        return out

    @staticmethod
    def convert_2_int(data):
        return map(lambda x: map(int, x), data)

    @staticmethod
    def remove_trailing_spaces(data):
        # assume two dimensional array
        return map(lambda x: map(str.strip, x), data)

    @staticmethod
    def remove_quotes_and_lower(data):
        data = map(lambda x: map(str.lower, x), data)
        return [[d[0].replace("'", "")] for d in data]

    @staticmethod
    def make_dic(featurenames):
        return {i: featurenames[i] for i in xrange(len(featurenames))}

    @staticmethod
    def make1D(data):
        if len(data[0]) == 1:
            return map(lambda x: x[0], data)
        else:
            print 'Cannot make it 1D!'
            sys.exit(1)

    @staticmethod
    def get_metafeature_to_orig_feature_mapping(metafeatures_array):
        metafeaturesDic = {}
        for i, elem in enumerate(metafeatures_array):
            metafeaturesDic[i] = [int(a) - 1 for a in elem]
        return metafeaturesDic

    @staticmethod
    def get_index_none_zeros(inarray):
        return [list(a.nonzero()[0]) for a in array(inarray)]
