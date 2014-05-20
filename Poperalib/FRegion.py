from __future__ import division
from __future__ import print_function
from numpy import *
from countreads import *
from multiprocessing import Pool
from countreads import *
import timeit
import sys

class KeyboardInterruptError(Exception):

    pass


class FRegion:

    def __init__(self, bamfile, nthreads, countchr=[]):

        self.bamfile = bamfile

        self.count_chr = countchr

        self.nthreads = nthreads

        self.__filte_region()

    def filte_region(self):

        bam_file = self.bamfile

        count_chr = self.count_chr

        nthreads = self.nthreads

        res = filter_region(bam_file, count_chr, nthreads)

        filted_region = res['filted_region']

        thresh_hold = res['thresh_hold']

        scare_std = res['region_std']

        scare_mean = res['region_mean']

        chr_total_reads = res['chr_total_reads']

        chrs_length = res['chrs_length']

        chrsfrcount = res['chrfrcount']

        filterreadscount = res['filterreadscount']

        totalreads = res['totalreads']

        chr_unique = res['chr_unique']

        adjreads = totalreads - filterreadscount

        self.filted_region = filted_region

        self.thresh_hold = thresh_hold

        self.region_std = scare_std

        self.region_mean = scare_mean

        self.chr_total_reads = chr_total_reads

        self.chrs_length = chrs_length

        self.chrsfcount = chrsfrcount

        self.totalreads = totalreads

        self.filterreadscount = filterreadscount

        self.chr_unique = chr_unique

        self.adjreads = adjreads

    __filte_region = filte_region


def filter_region(bamfile, count_chr, nthreads):

    samfile = pysam.Samfile(bamfile)

    chrs = list()

    filted_region = list()

    windowsize = 1000

    regionlist = list()

    refere_ncenumber = samfile.nreferences

    ref_lengths = samfile.lengths

    sam_ref = samfile.references

    chrs_length = dict()

    for chromosome in count_chr:

        chrs.append(chromosome)

    pars = list()

    chr_total_reads = dict()

    totalreads = 0

    for chromosome in count_chr:

        for i in range(refere_ncenumber):

            if sam_ref[i] == chromosome:

                chr_length = ref_lengths[i]

                chrs_length[chromosome] = chr_length

                par = dict()

                regionnow = chromosome+':'+str('1')+'-'+str(chr_length)

                par['region'] = regionnow

                par['chromosome'] = chromosome

                par['windowsize'] = windowsize

                par['bamfile'] = bamfile

                pars.append(par)

                chrregion = chromosome+':'+str('1')+'-'+str(chr_length)

                chrcount = dhsingleregioncounter(bamfile=bamfile, region=chrregion)

                chr_total_reads[chr] = chrcount

                totalreads = chrcount + totalreads

    pool = Pool(nthreads)

    try:

        windowcountlist = list()

        windowregionlist = list()

        starttime2 = timeit.timeit()

        chrswindow = pool.map(chrwindow_counter, pars)

        chr_unique = dict()

        for nowchrcount in chrswindow:

            nowchromosome = nowchrcount['chromosome']

            nowwindowcount = nowchrcount['windowcount']

            nowchruniq = nowchrcount['uniquesite']

            chr_unique[nowchromosome] = nowchruniq

            #for debug test =====
            # print (nowchromosome, nowchrcount['region'], "unique site", nowchruniq)
            #for debug test =====

            for nowscare in nowwindowcount:

                nowstart = nowscare * windowsize + 1

                nowend = (nowscare+1) * windowsize

                if nowend > chrs_length[nowchromosome]:

                    nowend = chrs_length[nowchromosome]

                nowregion = nowchromosome+":"+str(nowstart)+"-"+str(nowend)

                windowcountlist.append(nowwindowcount[nowscare])

                windowregionlist.append(nowregion)

        endtime2 = timeit.timeit()

        time2 = endtime2 - starttime2

        # print ("time2", time2)

        scare_mean = mean(windowcountlist)

        scare_std = std(windowcountlist)

        # print ("mean:", scare_mean, "std",scare_std)

        thresh_hold = scare_mean + 10 * scare_std

        chrsfrcount = 0

        filterreadscount = 0


        for i in range(0, len(windowcountlist)-1):

            if windowcountlist[i] >= thresh_hold:

                # print (windowregionlist[i]," reads count ", windowcountlist[i])

                filted_region.append(windowregionlist[i])

                filterreadscount = filterreadscount + windowcountlist[i]

        res = dict()

        res['filted_region'] = filted_region

        res['thresh_hold'] = thresh_hold

        res['region_std'] = scare_std

        res['region_mean'] = scare_mean

        res['chr_total_reads'] = chr_total_reads

        res['chrs_length'] = chrs_length

        res['chrfrcount'] = chrsfrcount

        res['filterreadscount'] = filterreadscount

        res['totalreads'] = totalreads

        res['chr_unique'] = chr_unique

        pool.close()

        return res

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


def chrwindow_counter(par):

    try:
        regionnow = par['region']

        chromosome = par['chromosome']

        windowsize = par['windowsize']

        bamfile = par['bamfile']

        windowcount = dhsinglewindowscarecounter(bamfile=bamfile, region=regionnow,
                                                 windowsize=windowsize)

        chrwindowcount = dict()

        chrwindowcount['windowcount'] = windowcount['readscount']

        chrwindowcount['chromosome'] = chromosome

        chrwindowcount['region'] = regionnow

        chrwindowcount['uniquesite'] = windowcount['uniquesite']

        return chrwindowcount

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)


if __name__ == "__main__":
    try:
        # main()
        import re
        from bgcount import *
        datafile = "Reb1_600mM_nova.bam"
        countchr=list()

        samfile = pysam.Samfile(datafile)

        sam_ref = samfile.references

        for i in sam_ref:

            countchr.append(i)

        excludechr = 'chrMito'

        if (excludechr):

            excludchr = excludechr.split(',')

            for chri in excludchr:

                if not chri in sam_ref:

                    print (chri, 'not in the %s file' % datafile)

                    print ("try to selcet exclude Chr from", end =" : ")

                    print (sam_ref, sep=",")

                    sys.exit(1)

                else:

                    j = 0

                    for n in countchr:

                        if chri == n:

                            del countchr[j]

                        j = j + 1

        k = 0

        digstart = re.compile('^\d')

        for m in countchr:

            if digstart.match(m):

                del countchr[k]

                print("skip chr:", m)

            k = k + 1

        a = FRegion(bamfile=datafile, nthreads=4, countchr=countchr)

        uniqreate = dhuniquerate(fregion=a)

        print (uniqreate)



    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)