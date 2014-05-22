from __future__ import division
from countreads import *
from kernel import *
import numpy as np


class KeyboardInterruptError(Exception):

    pass


def regionsmooth(bamfile, regionchromosome, regionstart, regionend, chr_length, kernelsize):

    try:

        # regionchromosome, sesite = region.split(':')

        # startsite, endsite = sesite.split('-')
        #
        # startsite = int(startsite)
        #
        # endsite = int(endsite)

        regionstart = regionstart - kernelsize*2

        regionend = regionend + kernelsize*2

        if regionstart < 1:

            regionstart = 1

        if regionend > chr_length:

            regionend = chr_length

        renewlength = regionstart - regionend + 1

        # smoothed_score = np.repeat(0, renewlength)

        # resizeregion = chromosome+":"+str(renewstart)+"-"+str(renewend)

        dhsinglecount = dhsinglereadscounter(bamfile=bamfile, regionchromosome=regionchromosome,
                                             regionstart=regionstart, regionend=regionend)

        kernelnow = smooth_kernel(kernelsize)

        readcount = list()

        kernel_score = list()

        for w in sorted(kernelnow):

            kernel_score.append(kernelnow[w])

        for n in range(regionstart, regionend+1):

            nowcount = 0

            if n in dhsinglecount:

                nowcount = dhsinglecount[n]

            readcount.append(nowcount)

        nowsmoothed = np.correlate(array(readcount), kernel_score, "same")

        outputscore = dict()

        outputscore['chromosome'] = regionchromosome

        outputscore['score'] = dict()

        # print (smoothed_score[0])

        for j in range(0, renewlength):

            nowsite = j + regionstart

            nowscore = nowsmoothed[j]

            if (regionstart<=nowsite<=regionend):

                outputscore['score'][nowsite] = nowscore

        return outputscore

    except KeyboardInterrupt:

        print ("You cancelled the program!")

        sys.exit(1)