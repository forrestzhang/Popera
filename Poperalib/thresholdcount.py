

from .kernelsmooth import *
from multiprocessing import Pool
from .FRegion import *

class KeyboardInterruptError(Exception):

    pass


def countthreshold(bamfile, count_chr, frgion, threshold, kernelsize=200, windowsize=100000, nthreads=4):

    pars = list()

    for chromosome in count_chr:

        for i in range(0, int(frgion.chrs_length[chromosome]/windowsize)+1):

            whether_filter = False

            nowstart = 1 + i * windowsize

            nowend = (i + 1) * windowsize

            for filtered_region in frgion.filted_region:

                filtered_region_chromosome, sesite = filtered_region.split(':')

                filtered_region_startsite, filtered_region_endsite = sesite.split('-')

                filtered_region_startsite = int(filtered_region_startsite)

                filtered_region_endsite = int(filtered_region_endsite)

                if filtered_region_chromosome == chromosome:

                    if nowstart<=filtered_region_startsite<=nowend:

                        whether_filter = True

                    if nowend<=filtered_region_endsite<=nowend:

                        whether_filter = True

            if whether_filter:

                print ("skip",chromosome,str(nowstart),str(nowend))

                continue

            par = dict()

            par['bamfile'] = bamfile

            par['kernelsize'] = kernelsize

            par['chr_length'] = frgion.chrs_length[chromosome]

            if (nowend > frgion.chrs_length[chromosome]):

                nowend = frgion.chrs_length[chromosome]

            nowregion = chromosome+":"+str(nowstart)+"-"+str(nowend)

            par['regionchromosome'] = str(chromosome)

            par['regionstart'] = nowstart

            par['regionend'] = nowend

            par['region'] = nowregion

            pars.append(par)

    pool=Pool(nthreads)

    try:

        score_threads = pool.map(countthresholdrunner, pars)

        score = dict()

        for nowscores in score_threads:

            for sc in nowscores:

                if sc in score:

                    score[sc] = score[sc] + nowscores[sc]

                else:

                    score[sc] = nowscores[sc]

        sumscore = 0

        totalbp = 0

        outtest = open("scoretest.txt",'w')

        for sc in sorted(score.keys()):
            print (str(sc)+"\t"+str(score[sc]),file=outtest)
            sumscore = score[sc] * sc + sumscore

            totalbp = score[sc] + totalbp
        outtest.close()

        meanscore = sumscore / totalbp

        sumsqure = 0

        for sc in score:

            sumsqure = sc*(score[sc] - meanscore)**2 + sumsqure

        std = (sumsqure/totalbp) ** (0.5)

        scorethreshold = meanscore + threshold * std

        print ("sumscore",sumscore,"totalbp",totalbp,"mean:", meanscore, "std:", std, "scorethreshold:",scorethreshold)

        pool.close()

        return scorethreshold

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


def countthresholdrunner(par):

    try:
        bamfile = par['bamfile']

        kernelsize = par['kernelsize']

        chr_length = par['chr_length']

        region = par['region']

        regionchromosome = par['regionchromosome']

        regionstart = par['regionstart']

        regionend = par['regionend']



        smoothscore = regionsmooth(bamfile, regionchromosome, regionstart, regionend, chr_length, kernelsize)

        chromosome, sesite = region.split(':')

        startsite, endsite = sesite.split('-')

        startsite = int(startsite)

        endsite = int(endsite)

        score=dict()

        for i in range(startsite, endsite):

            nowscore = 0

            if i in smoothscore['score']:

                nowscore = round(smoothscore['score'][i],3)

            if nowscore in score:

                score[nowscore] = score[nowscore] + 1

            else:

                score[nowscore] = 1

        return score

    except KeyboardInterrupt:

        raise KeyboardInterruptError()




if __name__ == "__main__":
    try:
        # main()
        import re
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

        a = FRegion(bamfile=datafile, countchr=countchr, nthreads=4)

        th = countthreshold(bamfile=datafile,count_chr=countchr,frgion=a, threshold=4)


        print ("thnow", th)



    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)