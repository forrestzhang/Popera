import sys
from decimal import Decimal, localcontext
from scipy.special import gamma, gammaincc
from scipy import math
import scipy.stats as stats

def bayesfactor(locallambda, peakscore):

    try:

        # bayesfactor = 2 * (math.log((gammaincc(peakscore-1, locallambda)*gamma(peakscore-1)), math.e) - (peakscore-1)*math.log(locallambda, math.e) + locallambda)
        #
        # a = (math.log(gammaincc(peakscore-1, locallambda), math.e) )
        # b = math.lgamma(peakscore-1)
        # c=(peakscore-1)*math.log(locallambda, math.e)
        # print (locallambda,peakscore,a,b,c)
        bayesfactor2 = 2 * (math.log(gammaincc(peakscore-1, locallambda), math.e)+math.lgamma(peakscore-1) - (peakscore-1)*math.log(locallambda, math.e) + locallambda)

        return bayesfactor2

    except Exception as e:

        print ('got exception in Jazzlib.sta.bayesfactor: %r,' % (e,))

        print (locallambda, peakscore)

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)



def poissonpvalue(x,mu):

    poissonpvalue = Decimal(1) - Decimal(stats.poisson.cdf(x, mu))

    return poissonpvalue



def fdr(pnow, plist, prank):
    #FDR=length(pvalue)*pvalue/rank(pvalue)

    rankofplist = prank

    lengthofplist = len(plist)

    for i in range(0,lengthofplist):

        if plist[i] == pnow:
            now_rank = rankofplist[i]
            fdr = lengthofplist*pnow/now_rank
            fdr = min(1,fdr)
            break

    return fdr
