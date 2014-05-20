from __future__ import division
from __future__ import print_function
from scipy.special import gamma, gammaincc
from scipy import math
import scipy.stats as stats

def bayesfactor(locallambda, peakscore):

    bayesfactor = 2 * (math.log((gammaincc(peakscore-1, locallambda)*gamma(peakscore-1)), math.e) - (peakscore-1)*math.log(locallambda, math.e) + locallambda)

    if 0 < bayesfactor == float("inf"):

        bayesfactor = 1419.565425786768

    return bayesfactor