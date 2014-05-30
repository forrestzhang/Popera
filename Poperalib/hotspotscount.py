from __future__ import division
from __future__ import print_function
from multiprocessing import Pool
from region import *
from numpy import *
from kernel import *
from countreads import *
import sys
from Hotspot import *
from bgcount import *


class KeyboardInterruptError(Exception):

    pass


def hotspotscount_nocontrol(bamfile, threshold, kernellength,
                    windowsize, nthreads, minlength,
                    fregion, samplename, countchr=[]):
    """
        hotspots[hid][type]

    """

    try:
        hotspots = list()

        uniqreate = dhuniquerate(fregion=fregion)

        cutoff = dhnoncontrol(uniqueratio=uniqreate, threshold=threshold, kernellength=kernellength, nthreads=nthreads)

        for chromosome in countchr:


            ultr = ultratio(uniqreate=uniqreate, fregion=fregion)

            chrlength = fregion.chrs_length[chromosome]

            chrhostspots = chrhotspotscount_nocontrol(bamfile, chromosome, ultr, cutoff,
                                                    kernellength, windowsize, nthreads, minlength, chrlength)

            # nowtype

            idnumber = 1

            for hotspotsnow in chrhostspots:

                start_site = hotspotsnow['start_site']

                end_site = hotspotsnow['end_site']

                idnumber = idnumber + 1

                hotspotid = chromosome + '.' + samplename +str(idnumber)

                hotspot = Hotspot(start=start_site, end=end_site, chromosome=chromosome,
                                  hotspotid=hotspotid)



                # print (hotspot.hotspotid,hotspot.chromosome) ####

                hotspots.append(hotspot)

        return hotspots

    except KeyboardInterrupt:

            raise KeyboardInterruptError()

def chrhotspotscount_nocontrol(bamfile, chromosome, ultratio, cutoff,
                   kernellength, windowsize, nthreads, minlength,chrlength):
    """

        jobtype = 'dh', 'nhsingle', 'nhpaired'

    """
    ercr = effectregion(chrlength=chrlength, windowsize=windowsize, bw= kernellength)

    pars = list()

    for scare in ercr:

        par = dict()

        ctstart = ercr[scare]['ctstart']

        ctend = ercr[scare]['ctend']

        par['ultratio'] = ultratio

        par['chromosome'] = chromosome

        par['regionstart']= ctstart

        par['regionend']= ctend

        par['kernellength']= kernellength

        par['bamfile'] = bamfile

        par['scare'] = scare

        par['cutoff'] = cutoff

        pars.append(par)

    pool=Pool(nthreads)

    chrhotspots = dict()

    try:

        hotsports_point = pool.map(hostspotspointcounter_nocontrol, pars)

        for hotspots_scare in hotsports_point:

            for nowsite in hotspots_scare:

                if not nowsite in chrhotspots:

                    if (nowsite >= 0 and nowsite<chrlength):

                        chrhotspots[nowsite] = 1

        pool.close()

        del hotsports_point

        hotspots_list=chrhotspots.keys()

        chrhotsopts_region = continueregion(points=hotspots_list, minlength=minlength)

        return chrhotsopts_region

    except KeyboardInterrupt:

        pool.terminate()

        print ("You cancelled the program!")

        sys.exit(1)

    except Exception, e:

        print ('got exception: %r, terminating the pool' % (e,))

        pool.terminate()

        print ('pool is terminated')

    finally:
    #     print ('joining pool processes')
        pool.join()
        # print ('join complete')


def hostspotspointcounter_nocontrol(par):

    try:

        #ctstart = ercr[scare]['ctstart']
        #ctend = ercr[scare]['ctend']

        ultratio = par['ultratio']

        chromosome = par['chromosome']

        regionstart = par['regionstart']

        regionend = par['regionend']

        kernellength = par['kernellength']

        bamfile = par['bamfile']

        scare = par['scare']

        cutoff = par['cutoff']

        ###get par

        kernel = smooth_kernel(kernellength)

        samfile = pysam.Samfile(bamfile)

        # countregion = chromosome + ':' + str(regionstart) + '-' + str(regionend)

        hotsopts_scare = list()

        reads_score = dict()

        reads_score = dhsinglereadsnormailzed(bamfile=bamfile, regionchromosome=chromosome, regionstart=regionstart,
                                              regionend=regionend, ultratio=ultratio)

        kernel_score = list()

        for i in sorted(kernel):

            kernel_score.append(kernel[i])

        readsscore_list = list()
        # print (countregion,regionstart, regionend) ###
        for n in range(regionstart, regionend):

            nowscore = 0

            if reads_score.has_key(n):

                nowscore = reads_score[n]

            readsscore_list.append(nowscore)

        smoothed_result = correlate(array(readsscore_list), kernel_score, "same")

        for i in range(0, smoothed_result.size):

            if smoothed_result[i] > cutoff:

                now_site = i + regionstart

                hotsopts_scare.append(now_site)

        return hotsopts_scare

    except KeyboardInterrupt:

        raise KeyboardInterruptError()