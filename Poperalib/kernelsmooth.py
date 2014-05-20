from __future__ import division
from countreads import *
from kernel import *
import numpy as np


class KeyboardInterruptError(Exception):

    pass


def regionsmooth(bamfile, region, chr_length, kernelsize):

    chromosome, sesite = region.split(':')

    startsite, endsite = sesite.split('-')

    startsite = int(startsite)

    endsite = int(endsite)

    renewstart = startsite - kernelsize*2

    renewend = endsite + kernelsize*2

    if renewstart < 1:

        renewstart = 1

    if renewend > chr_length:

        renewend = chr_length

    renewlength = renewend - renewstart + 1

    # smoothed_score = np.repeat(0, renewlength)

    resizeregion = chromosome+":"+str(renewstart)+"-"+str(renewend)

    dhsinglecount = dhsinglereadscounter(bamfile=bamfile, region=resizeregion)

    kernelnow = smooth_kernel(kernelsize)

    readcount = list()

    kernel_score = list()

    for w in sorted(kernelnow):

        kernel_score.append(kernelnow[w])

    for n in range(renewstart, renewend+1):

        nowcount = 0

        if n in dhsinglecount:

            nowcount = dhsinglecount[n]

        readcount.append(nowcount)

    nowsmoothed = np.correlate(array(readcount), kernel_score, "same")

    outputscore = dict()

    outputscore['chromosome'] = chromosome

    outputscore['score'] = dict()

    # print (smoothed_score[0])

    for j in range(0, renewlength):

        nowsite = j + renewstart

        nowscore = nowsmoothed[j]

        if (startsite<=nowsite<=endsite):

            outputscore['score'][nowsite] = nowscore

    return outputscore
