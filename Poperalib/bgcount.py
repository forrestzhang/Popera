

import pysam
from numpy import *
from multiprocessing import Pool
import random as rnd
from .kernel import *
from .countreads import *
import sys
from .FRegion import *

class KeyboardInterruptError(Exception):

    pass


def dhnoncontrol(uniqueratio, threshold, kernellength, nthreads=4):

    bgscore = sim_replicate_nthreads(run_times=200, uniqueratio=uniqueratio,
                                     nthreads=nthreads, kernellength=kernellength, threshold=threshold)

    cutoff = bgscore['mean'] + bgscore['std'] * threshold

    return cutoff


def sim_replicate_nthreads(run_times=1000, uniqueratio=1, kernellength = 600, threshold = 4, nthreads = 2):
    # randomthresh = list()

    pars = list()

    for i in range(0,run_times):

        par=dict()

        par['uniqueratio'] = uniqueratio

        par['kernellength'] = kernellength

        par['threshold'] = threshold

        pars.append(par)

    pool=Pool(nthreads)

    outscore = dict()

    try:
        randomthresh = pool.map(sim_bg_thread_worker, pars)

        summean = 0.0

        sumstd = 0.0

        for randscore in randomthresh:

            randmean = randscore['rand_mean']

            randstd = randscore['rand_std']
            # print (randmean, randstd)
            summean = summean + randmean

            sumstd = sumstd + randstd

        mean_of_mean = summean/run_times

        mean_of_std = sumstd/run_times
        # print ('mean_of_mean',mean_of_mean, 'mean_of_std',mean_of_std)

        outscore['mean'] = mean_of_mean

        outscore['std'] = mean_of_std
        #return (mean_of_mean, mean_of_std)

        pool.close()

        return outscore

    except KeyboardInterrupt:

        pool.terminate()

        print ("You cancelled the program!")

        sys.exit(1)

    except Exception as e:

        print ('got exception: %r, terminating the pool' % (e,))

        pool.terminate()

        print ('pool is terminated')

    finally:
        # print ('joining pool processes')
        pool.join()
    # print ('join complete')
    # pool.join()
    # pool.close()


def sim_bg_thread_worker(par):

    try:

        uniqueratio=par['uniqueratio']

        kernellength = par['kernellength']

        threshold = par['threshold']

        kernel = smooth_kernel(length=kernellength)

        sim_genome_size = int(1e5)

        total_reads = int(sim_genome_size * uniqueratio)

        region_site = list(range(0,sim_genome_size))

        sim_uniqsite = rnd.sample(region_site, total_reads)

        rand_reads_count = list()

        for i in range(0,sim_genome_size):

            rand_reads_count.append(0)

        kernel_score = list()

        for i in sorted(kernel):
            kernel_score.append(kernel[i])



        kdesmooth_result = dict()

        for i in range(0,total_reads):

            rand_number = int(rnd.uniform(0,total_reads))

            rand_reads = sim_uniqsite[rand_number]

            rand_reads_count[rand_reads] = rand_reads_count[rand_reads] + 1.0

        smoothed_result = correlate(array(rand_reads_count), kernel_score)

        scores = list()

        rand_mean = smoothed_result.mean()

        rand_std = smoothed_result.std()

        total_sum = smoothed_result.sum()

        rand_threshhold = rand_mean + threshold * rand_std

        higher_count = 0

        for now_site in kdesmooth_result:

            if kdesmooth_result[now_site] > rand_threshhold:

                higher_count = higher_count + 1

        # print (total_sum, rand_mean, rand_std, rand_threshhold, higher_count, total_reads)

        randscore = dict()

        randscore['rand_mean'] = rand_mean

        randscore['rand_std'] = rand_std

        return randscore

    except KeyboardInterrupt:

        raise KeyboardInterruptError()


def dhuniquerate(fregion):

    uniqsite = 0

    genomelength = 0

    for chromosome in fregion.chrs_length:

        genomelength = genomelength + int(fregion.chrs_length[chromosome])

    for chromosome in fregion.chr_unique:

        uniqsite = uniqsite + int(fregion.chr_unique[chromosome])

    # print (fregion.chr_unique)
    #
    # print (fregion.chrs_length)

    uniqreate = uniqsite/genomelength

    return uniqreate


def ultratio(fregion, uniqreate):
    """
    ultraio = chrlength * uniqueratio / chr_total_reads
    """
    genomelength = 0

    adjreads = fregion.adjreads

    for chromosome in fregion.chrs_length:

        genomelength = genomelength + int(fregion.chrs_length[chromosome])

    ultratio = genomelength * uniqreate / adjreads

    return ultratio

def get_bpc(bamfile,  hotspots, filted_region, nthreads, maxinsert = 100000):

    #bpc average readscount per basepare

    # total_reads = 0
    #
    total_length = 0
    #
    # samfile = pysam.Samfile(bamfile)

    pars = list()

    for hotspot_now in hotspots:

        par = dict()

        par['bamfile'] = bamfile

        par['hotspot'] = hotspot_now

        par['filted_region'] = filted_region

        par['maxinsert'] = maxinsert

        pars.append(par)

        total_length = hotspot_now.end - hotspot_now.start + 1 + total_length


    # print ("total length %s" % total_length)

    pool = Pool(nthreads)

    try:

        reads_count = pool.map(bpc_runner, pars)

        total_reads = 0.0

        for count_now in reads_count:

            total_reads = total_reads + count_now

        bpc = (total_reads+0.0)/total_length

        pool.close()

        return bpc

    except KeyboardInterrupt:

        pool.terminate()

        print ("You cancelled the program!")

        sys.exit(1)

    except Exception as e:

        print ('got exception: %r, terminating the pool' % (e,))

        pool.terminate()

        print ('pool is terminated')

    finally:
        # print ('joining pool processes')
        pool.join()
        # print ('join complete')


def bpc_runner(par):

    try:

        bamfile = par['bamfile']

        hotspot = par['hotspot']

        filted_region = par['filted_region']

        maxinsert = par['maxinsert']

        start_site = hotspot.start

        end_site = hotspot.end

        whether_in_fr = 0

        chromosome = hotspot.chromosome

        #hotspotregio = chromosome + ':' + str(start_site) + '-' + str(end_site)

        hotspotreads = 0

        for i in range(start_site, end_site + 1):

            parentscare = int(i/100)

            if chromosome in filted_region:

                if parentscare in filted_region[chromosome]:

                    whether_in_fr = 1

        if whether_in_fr == 0:

            readscount = dict()

            readscount = dhsinglereadscounter(bamfile = bamfile, regionchromosome = chromosome, 
                                              regionstart = start_site, regionend = end_site)

            for i in readscount:

                hotspotreads = hotspotreads+readscount[i]

        # print (hotspotreads)
        return hotspotreads

    except KeyboardInterrupt:

        raise KeyboardInterruptError()