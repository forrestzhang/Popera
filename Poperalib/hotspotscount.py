

from multiprocessing import Pool
from .region import *
from numpy import *
from .kernel import *
from .countreads import *
import sys
from .Hotspot import *
from .bgcount import *
# from cEM_zip import *
from .sta import *

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

        hotspots_list=list(chrhotspots.keys())

        chrhotsopts_region = continueregion(points=hotspots_list, minlength=minlength)

        return chrhotsopts_region

    except KeyboardInterrupt:

        pool.terminate()

        print ("You cancelled the program!")

        sys.exit(1)

    except Exception as e:

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

            if n in reads_score:

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


def hotspots_bayes(bayesfactorthreshold,nthreads, hotspots, bamfile, fregion, exsize):

    pool = Pool(nthreads)

    try:

        pars = list()

        for hotspot in hotspots:

            par = dict()

            par['hotspt'] = hotspot

            par['datafile'] = bamfile

            par['exsize'] = exsize

            par['fregion'] = fregion

            par['bayesfactorthreashold'] = bayesfactorthreshold

            pars.append(par)

        bayeshotspots = pool.map(hotspots_bayes_worker, pars)

        newhotspots = list()

        for hotspotnow in bayeshotspots:

            newhotspots.append(hotspotnow)

        pool.close()

        return newhotspots

    except KeyboardInterrupt:

        pool.terminate()

        print ("You cancelled the program!")

        sys.exit(1)

    except Exception as e:

        print ('got exception in Poperalib.hotspotscount.hotspots_bayes: %r, terminating the pool' % (e,))

        pool.terminate()

        print ('pool is terminated')

    finally:

        pool.join()


# def hotspots_bayes_worker(par):
#
#     try:
#
#         hotspot = par['hotspt']
#
#         bamfile = par['datafile']
#
#         exsize = par['exsize']
#
#         fregion = par['fregion']
#
#         bayesfactorthreshold = par['bayesfactorthreashold']
#
#         start = hotspot.start
#
#         end = hotspot.end
#
#         chromosome = hotspot.chromosome
#
#         chrlength = fregion.chrs_length[chromosome]
#
#         regionstart = start - 5100
#
#         regionend = end + 5100
#
#         if regionstart < 1:
#
#             regionstart = 1
#
#         if regionend > chrlength:
#
#             regionend = chrlength
#
#         enrichedsite = dict()
#
#         bayesfactorscore = dict()
#
#         inputwindow5k = list()
#
#         inputwindow10k = list()
#
#         datacount = dhsinglereadsexcount(bamfile=bamfile, regionchromosome=chromosome, regionstart=regionstart,
#                                          regionend=regionend, exsize=exsize)
#
#         #print (datacount)
#
#         for sitenow in range(start-5000, end+5000):
#
#             nowcount = 0
#
#             if sitenow < 0:
#
#                 continue
#
#             if sitenow > chrlength:
#
#                 continue
#
#             if sitenow in datacount:
#
#                 nowcount = datacount[sitenow]
#
#             inputwindow10k.append(nowcount)
#
#         for sitenow in range(start-2500,end+2500):
#
#             nowcount = 0
#
#             if sitenow < 0:
#
#                 continue
#
#             if sitenow > chrlength:
#
#                 continue
#
#             if sitenow in datacount:
#
#                 nowcount = datacount[sitenow]
#
#             inputwindow5k.append(nowcount)
#
#
#         (window5klhat, window5kphat) = cEM_zip(inputwindow5k)
#
#         (window10klhat, window10kphat) = cEM_zip(inputwindow10k)
#
#         maxlhat = max(window5klhat, window10klhat)
#
#         maxbayes = 0
#
#         maxsite = 0
#
#         for wsite in range(start-1, end+1):
#
#             if wsite in datacount:
#
#                 nowcount = int(datacount[wsite])
#
#                 if nowcount < maxlhat:
#
#                     continue
#
#                 elif nowcount < 2:
#
#                     continue
#
#                 else:
#
#                     nowbayesfactor = bayesfactor(locallambda=maxlhat, peakscore=nowcount)
#
#
#                     if nowbayesfactor > maxbayes:
#
#                         maxbayes = nowbayesfactor
#
#                         maxsite = wsite
#
#         hotspot.bayescore = maxbayes
#
#         # print (hotspot.hotspotid,hotspot.bayescore,maxbayes)
#
#         return hotspot
#
#
#     except Exception as e:
#
#         print ('got exception in Poperalib.hotspotscount.hotspots_bayes_worker: %r, terminating the pool' % (e,))
#
#         print (par['hotspot'].chromosome, par['hotspot'].start,par['hotspot'].end)
#
#
#     except KeyboardInterrupt:
#
#         print ("You cancelled the program!")
#
#         sys.exit(1)