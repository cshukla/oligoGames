#!/usr/bin/env python
from optparse import OptionParser
import numpy as np
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression, Ridge, RidgeCV
from sklearn.cross_decomposition import PLSRegression
from sklearn.svm import SVR
from sklearn.gaussian_process import GaussianProcess
from sklearn.cross_validation import KFold
import copy, sys
import dna

################################################################################
# sequence_regression.py
#
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fasta> <scores>'
    parser = OptionParser(usage)
    parser.add_option('-a', dest='canonical_kmers', default=False, action='store_true', help='Count canonical k-mers [Default: %default]')
    parser.add_option('--alpha', dest='alpha', default=None, type='float', help='Regularization alpha parameter. Will choose via CV if not specified [Default: %default]')
    parser.add_option('-c', dest='cv_folds', default=0, type='int', help='Cross-validate with this many folds [Default: %default]')
    parser.add_option('--epsilon', dest='epsilon', default=None, type='float', help='Regularization epsilon parameter. Will choose via CV if not specified [Default: %default]')
    parser.add_option('-g', dest='gaps', default=0, type='int', help='Gaps in k-mers string kernel [Default: %default]')
    parser.add_option('-k', dest='k', default=4, type='int', help='K-mer size for string kernel [Default: %default]')
    parser.add_option('-l', dest='length', default=False, action='store_true', help='Add log2 sequence length as an attribute [Default: %default]')
    parser.add_option('-m', dest='method', default='ols', help='Regression method [Default: %default]')
    parser.add_option('-o', dest='output_file', default='seq_regr.txt', help='Output file [Default: %default]')
    parser.add_option('-w', dest='whiten', default=False, action='store_true', help='Whiten the sequence scores [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide fasta file and scores file')
    else:
        fasta_file = args[0]
        scores_file = args[1]

    ##################################################
    # convert sequences to feature representations
    ##################################################
    seq_vectors = fasta_string_kernel(fasta_file, options.k, options.gaps, options.canonical_kmers)

    if options.length:
        add_length_feature(seq_vectors, fasta_file)

    ##################################################
    # read scores
    ##################################################
    seq_scores = {}

    scores_in = open(scores_file)
    
    try: 
        line = scores_in.readline()
        a = line.split()
        seq_scores[a[0]] = float(a[1])
    except:
        # possible header line
        pass
        
    for line in scores_in:
        a = line.split()
        seq_scores[a[0]] = float(a[1])


    ##################################################
    # make scikit-learn data structures
    ##################################################
    # shitty method filling in the dense matrix
    kmers = set()
    for kmer_vec in seq_vectors.values():
        kmers |= set(kmer_vec.keys())

    kmers_sort = sorted(kmers)

    seq_headers = sorted(seq_vectors.keys())
    
    X = np.array([[seq_vectors[header].get(kmer,0) for kmer in kmers_sort] for header in seq_headers])
    y = np.array([seq_scores[header] for header in seq_headers])

    if options.whiten:
        y = preprocessing.scale(y)

    ##################################################
    # decide method
    ##################################################
    if options.method.lower() == 'ols':
        model = LinearRegression()

    elif options.method.lower() == 'pls':
        model = PLSRegression(n_components=2)

    elif options.method.lower() == 'ridge':
        if options.alpha:
            # model = Ridge(alpha=options.alpha)
            model = RidgeCV(alphas=[options.alpha], store_cv_values=True)
        else:
            #model = RidgeCV(alphas=[0.0001, 0.0002, 0.0004, 0.0008, .0016, 0.0032, 0.0064, .0128], store_cv_values=True)
            model = RidgeCV(alphas=[0.0004, 0.0008, 0.0016, 0.0032], store_cv_values=True)

    elif options.method.lower() == 'svm':
        if options.alpha:
            svm_c = len(y) / options.alpha
        else:
            svm_c = 100
        if options.epsilon:
            svm_eps = options.epsilon
        else:
            svm_eps = 0.5

        model = SVR(kernel='linear', degree=3, C=svm_c, epsilon=svm_eps)

    elif options.method.lower() == 'gp':
        model = GaussianProcess()

    else:
        print >> sys.stderr, 'Method not recognized.'
        exit(1)


    ##################################################
    # learn model
    ##################################################
    model.fit(X, y)

    ss_tot = sum(np.square(y - np.mean(y)))

    if options.method.lower() == 'ridge':
        for i in range(len(model.alphas)):
            score_cv = (1.0 - sum(model.cv_values_[:,i])/ss_tot)
            print >> sys.stderr, 'RidgeCV alpha=%.5f score=%f' % (model.alphas[i], score_cv)

    ##################################################
    # cross-validate
    ##################################################
    if options.cv_folds > 0:
        scores = []
        ss_reg = 0

        if options.method.lower() == 'ridge':
            model_cv = Ridge(alpha=model.alpha_)
        else:
            model_cv = copy.copy(model)
        
        kf = KFold(len(y), n_folds=options.cv_folds, shuffle=True)
        for train, test in kf:
            X_train, X_test, y_train, y_test = X[train], X[test], y[train], y[test]

            # learn on train
            model_cv.fit(X[train], y[train])

            # score on test
            scores.append(model_cv.score(X_test, y_test))

            ss_reg += sum(np.square(y_test - model_cv.predict(X_test)))

        score_cv = 1 - ss_reg / ss_tot

            
    ##################################################
    # output model information
    ##################################################
    model_out = open(options.output_file, 'w')

    print >> model_out, 'Score\t%.3f' % model.score(X, y)
    if options.cv_folds > 0:
        print >> model_out, 'ScoreCV\t%.3f' % score_cv
        if options.method.lower() == 'ridge' and options.alpha:
            score_cv = (1.0 - sum(model.cv_values_)/ss_tot)
            print >> model_out, 'ScoreCV\t%.3f' % score_cv

    for i in range(len(kmers_sort)):
        if options.method.lower() == 'pls':
            coef_i = model.coefs[i]
        else:
            coef_i = model.coef_[i]

        print >> model_out, '%s\t%f' % (kmers_sort[i], coef_i)

    model_out.close()
    


################################################################################
# add_length_feature
#
# Add log2 sequence length as a feature.
################################################################################
def add_length_feature(seq_vectors, fasta_file):
    seq_lengths = {}
    for line in open(fasta_file):
        if line[0] == '>':
            header = line[1:].rstrip()
            seq_lengths[header] = 0
        else:
            seq_lengths[header] += len(line.rstrip())
    
    for header in seq_lengths:
        seq_vectors[header]['length'] = np.log2(seq_lengths[header])


################################################################################
# fasta_string_kernel
#
# Compute a string kernel profile for each sequence in the fasta file.
################################################################################
def fasta_string_kernel(fasta_file, k, gaps, canonical):
    seq_vectors = {}
    seq = ''
    for line in open(fasta_file):
        if line[0] == '>':
            if seq:
                if gaps == 0:
                    seq_vectors[header] = kmer_kernel(seq, k, canonical)
                elif gaps == 1:
                    seq_vectors[header] = kmer_mismatch1_kernel(seq, k)
                else:
                    print >> sys.stderr, 'Gaps >1 not implemented'
                    exit(1)

            header = line[1:].rstrip()
            seq = ''
        else:
            seq += line.rstrip()

    if gaps == 0:
        seq_vectors[header] = kmer_kernel(seq, k, canonical)
    elif gaps == 1:
        seq_vectors[header] = kmer_mismatch1_kernel(seq, k)
    else:
        print >> sys.stderr, 'Gaps >1 not implemented'
        exit(1)

    return seq_vectors


################################################################################
# kmer_kernel
#
# Compute k-mer profile of seq into a dict.
################################################################################
def kmer_kernel(seq, k, canonical=True):
    kmer_counts = {}

    if canonical:
        seq_rc = dna.rc(seq)

    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]
        kmer_counts[kmer] = kmer_counts.get(kmer,0) + 1

        if canonical:
            kmer_rc = seq_rc[i:i+k]
            kmer_counts[kmer_rc] = kmer_counts.get(kmer_rc,0) + 1

    if canonical:
        kmer_counts = dna.canonical_kmers(kmer_counts)

    # normalize
    # kmer_sum = float(sum(kmer_counts.values()))
    kmer_sum = float(sum(np.square(kmer_counts.values())))

    vec = {}
    for kmer in kmer_counts:
        vec[kmer] = kmer_counts[kmer] / kmer_sum

    return vec


################################################################################
# kmer_mismatch1_kernel
#
# Compute 1 mismatch k-mer profile of seq into a dict.
################################################################################
def kmer_mismatch1_kernel(seq, k):
    vec = {}

    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]
        for j in range(k):
            kmer1 = kmer[:j] + '.' + kmer[j+1:]
            vec[kmer1] = vec.get(kmer1,0) + 1

    # normalize
    kmer_sum = float(sum(np.square(vec.values())))
    for kmer in vec:
        vec[kmer] /= kmer_sum

    return vec


################################################################################
# kmer_mismatch_kernel
#
# Compute g mismatch k-mer profile of seq into a dict.
################################################################################
def kmer_mismatch_kernel(seq, k, g):
    vec = {}

    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]

        gap_set = range(g)
        while True:
            gkmer = ''.join([kmer[j] if j in gap_set else '.' for j in range(k)])
            vec[gkmer] = vec.get(gkmer,0) + 1

            u = g-1
            while u >= 0:
                if gap_set[u] == k-1 or gap_set[u] == gap_set[u+1] - 1:
                    u -= 1
                else:
                    gap_set[u] += 1
                    # change the suffix

    # normalize
    kmer_sum = float(sum(np.square(kmer_counts.values())))

    for kmer in vec:
        vec[kmer] /= kmer_sum

    return vec



################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
