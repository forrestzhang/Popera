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

        startsite = regionstart

        endsite = regionend

        regionstart = regionstart - kernelsize*2

        regionend = regionend + kernelsize*2

        if regionstart < 1:

            regionstart = 1

        if regionend > chr_length:

            regionend = chr_length

        renewlength = regionend - regionstart + 1

        # print (regionstart, regionend)

        smoothed_score = np.repeat(0, renewlength)

        # resizeregion = chromosome+":"+str(renewstart)+"-"+str(renewend)

        dhsinglecount = dhsinglereadscounter(bamfile=bamfile, regionchromosome=str(regionchromosome),
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

            if (startsite<=nowsite<=endsite):

                outputscore['score'][nowsite] = nowscore

                # print (regionchromosome,nowsite,nowscore)

        return outputscore

    except KeyboardInterrupt:

        print ("You cancelled the program!")

        sys.exit(1)