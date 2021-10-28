

import io
from .Sampleinfor import *
from .countreads import *
from .Hotspot import *
from .bgcount import *
from multiprocessing import Pool
from .kernelsmooth import *
from .kernel import *
import numpy as np
import pyBigWig


class KeyboardInterruptError(Exception):

    pass


def normalizeratio(sampleinfors):

    normailziedratio = dict()

    adjreads = list()

    for sampleinfor in sampleinfors:

        adjreads.append(sampleinfor.fregion.adjreads)

    # minreadscount = min(adjreads)

    for sampleinfor in sampleinfors:

        normailziedratio[sampleinfor.samplename] = sampleinfor.fregion.adjreads/sampleinfor.fregion.countgenomelength

    return normailziedratio


def hotspotswriter(samplename, hotspots, bayesfactorthreshold=0):

    bedfilename =samplename+ '_' + 'hotspots' + ".bed"

    # open_bed = io.FileIO(bedfilename, 'w')
    open_bed = open(bedfilename, 'w')
    if bayesfactorthreshold == 0:

        for hotspot in hotspots:

            #bedlist = [str(hotspot.chromosome), str(hotspot.start), str(hotspot.end), hotspot.hotspotid]
            
            bedlist = [str(hotspot.chromosome), str(hotspot.start), str(hotspot.end), hotspot.peakid, str(hotspot.pvalue),'.']

            linker = "\t"

            outstring = linker.join(bedlist)

            # open_bed.write(outstring)
            print(outstring, file=open_bed)
    else:

        for hotspot in hotspots:

            #bedlist = [str(hotspot.chromosome), str(hotspot.start), str(hotspot.end), hotspot.hotspotid]
            if hotspot.bayescore >= bayesfactorthreshold:

                bedlist = [str(hotspot.chromosome), str(hotspot.start), str(hotspot.end), hotspot.peakid, str(hotspot.bayescore),float(hotspot.pvalue)]

                linker = "\t"

                outstring = linker.join(bedlist) + "\n"

                print(outstring, file=open_bed)

    open_bed.close()

def narrowpeakwriter(samplename, hotspots, bayesfactorthreshold=0):

    bedfilename =samplename+ '_' + 'hotspots' + ".narrowPeak"

    # open_bed = io.FileIO(bedfilename, 'w')
    open_bed = open(bedfilename, 'w')
    if bayesfactorthreshold == 0:

        for hotspot in hotspots:

            #bedlist = [str(hotspot.chromosome), str(hotspot.start), str(hotspot.end), hotspot.hotspotid]
            
            bedlist = [str(hotspot.chromosome), str(hotspot.start), str(hotspot.end), hotspot.peakid, str(hotspot.pvalue),'.','-1',str(hotspot.pvalue),'-1','-1']

            linker = "\t"

            outstring = linker.join(bedlist)

            # open_bed.write(outstring)
            print(outstring, file=open_bed)
    open_bed.close()

def mergedhotsportscountwrite(sampleinfors, mergedhotspots, nthreads, outfilename="merged_hotsopts_count.txt"):

    try:

        outfile = open(outfilename,'w')

        headerstring = "location"

        for sampleinfor in sampleinfors:

            headerstring = headerstring + "\t" + sampleinfor.samplename

        print(headerstring, file=outfile)

        normalizedratio = normalizeratio(sampleinfors=sampleinfors)

        pars = list()

        for hotspot in mergedhotspots:

            par = dict()

            par['hotspot'] = hotspot

            par['sampleinfors'] = sampleinfors

            par['normalizedratio'] = normalizedratio

            pars.append(par)

        pool = Pool(nthreads)



        countinthreads = pool.map(hotspotscounter, pars)

        for countstr in countinthreads:

            print (countstr, file=outfile)

        pool.close()

        outfile.close()

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

def singlehotsportcountwrite(sampleinfor,hotspots,nthreads):

    try:

        samplename = sampleinfor.samplename

        outfile = open('%s_hotspot_count.txt'%samplename,'w')

        print('location',samplename,sep='\t',file=outfile)

        normalizedratio = normalizeratio(sampleinfors=[sampleinfor])

        pars = []

        for hotspot in hotspots:

            par = dict()

            par['hotspot'] = hotspot

            par['sampleinfors'] = [sampleinfor]

            par['normalizedratio'] = normalizedratio

            pars.append(par)

        pool = Pool(nthreads)

        countinthreads = pool.map(hotspotscounter, pars)

        for countstr in countinthreads:

            print (countstr, file=outfile)

        pool.close()

        outfile.close()
    
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
    
def hotspotscounter(par):

    try:

        hotspot = par['hotspot']

        sampleinfors = par['sampleinfors']

        outstr = hotspot.region

        normalizedratio = par['normalizedratio']

        for sampleinfor in sampleinfors:

            nowreads = dhsingleregioncounter(bamfile=sampleinfor.datafile, regionchromosome=hotspot.chromosome,
                                             regionstart=hotspot.start, regionend=hotspot.end)

            normailziedcount = nowreads/(hotspot.end-hotspot.start+1)/normalizedratio[sampleinfor.samplename]

            outstr = outstr + "\t" + str(normailziedcount)

        return outstr

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)


def wigwritte(sampleinfors, kernellength, nthreads):

    try:

        normalizedratio = normalizeratio(sampleinfors=sampleinfors)

        pars = list()

        for sampleinfor in sampleinfors:

            uniqreate = dhuniquerate(fregion=sampleinfor.fregion)

            maxscore = dhnoncontrol(uniqueratio=uniqreate, threshold=50, kernellength=kernellength, nthreads=nthreads)

            samplenormalizedratio = normalizedratio[sampleinfor.samplename]

            par = dict()

            # par['uniqreate'] = uniqreate

            par['maxscore'] = maxscore

            par['sampleinfor'] = sampleinfor

            par['samplenormalizedratio'] = samplenormalizedratio

            par['kernellength'] = kernellength

            par['maxscore'] = maxscore/samplenormalizedratio

            pars.append(par)

        pool = Pool(nthreads)

        pool.map(wigwritter, pars)

        pool.close()

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


def wigwritter(par):

    try:
        sampleinfor = par['sampleinfor']

        kernellength = par['kernellength']

        samplenormalizedratio = par['samplenormalizedratio']

        kernellength = par['kernellength']

        maxscore = par['maxscore']

        bamfile = sampleinfor.datafile

        wigfilename = sampleinfor.samplename+".wig"

        wigio = open(wigfilename, 'w')

        kernel = smooth_kernel(kernellength)

        kernel_score = list()

        for w in sorted(kernel):

            kernel_score.append(kernel[w])

        print('track type=wiggle_0 name="',sampleinfor.samplename,'" description="', sampleinfor.samplename,'"',sep='',file=wigio)

        for chromosome in sampleinfor.fregion.count_chr:

            chromosome = str(chromosome)



            print('variableStep	chrom=', chromosome, sep='', file=wigio)

            # chrregion = chromosome+":"+str(1)+"-"+str(sampleinfor.fregion.chrs_length[chromosome])
            for scare in range(0, int(sampleinfor.fregion.chrs_length[chromosome]/1000000)+1):

                startsite = scare * 1000000 + 1

                endsite = (scare + 1) * 1000000

                if endsite > sampleinfor.fregion.chrs_length[chromosome]:

                    endsite = sampleinfor.fregion.chrs_length[chromosome]

                # print (chromosome,startsite,endsite)

                # regionnow = chromosome+":" + str(startsite) + "-" + str(endsite)

                smoothedscore = regionsmooth(bamfile=bamfile, regionchromosome=str(chromosome),
                                             regionstart=startsite, regionend=endsite,
                                             chr_length=sampleinfor.fregion.chrs_length[chromosome],
                                             kernelsize=kernellength)

                # print (smoothedscore['score'])

                for site in sorted(smoothedscore['score'].keys()):

                    # print (site, smoothedscore['score'][site])

                    score = smoothedscore['score'][site]/samplenormalizedratio

                    if score > maxscore:

                        score = maxscore

                    score = round(score, 3)

                    print (str(site)+"\t"+str(score), file=wigio)

        wigio.close()

    except KeyboardInterrupt:

        print ("You cancelled the program!")

        sys.exit(1)

def bigwigwritte(sampleinfors, kernellength, nthreads):

    try:

        normalizedratio = normalizeratio(sampleinfors=sampleinfors)

        pars = list()

        for sampleinfor in sampleinfors:

            uniqreate = dhuniquerate(fregion=sampleinfor.fregion)

            maxscore = dhnoncontrol(uniqueratio=uniqreate, threshold=50, kernellength=kernellength, nthreads=nthreads)

            samplenormalizedratio = normalizedratio[sampleinfor.samplename]

            par = dict()

            # par['uniqreate'] = uniqreate

            par['maxscore'] = maxscore

            par['sampleinfor'] = sampleinfor

            par['samplenormalizedratio'] = samplenormalizedratio

            par['kernellength'] = kernellength

            par['maxscore'] = maxscore/samplenormalizedratio

            pars.append(par)

        pool = Pool(nthreads)

        pool.map(bigwigwritter, pars)

        pool.close()

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


def bigwigwritter(par):

    try:
        sampleinfor = par['sampleinfor']

        kernellength = par['kernellength']

        samplenormalizedratio = par['samplenormalizedratio']

        kernellength = par['kernellength']

        maxscore = par['maxscore']

        bamfile = sampleinfor.datafile

        bwfilename = sampleinfor.samplename+".bw"

        #wigio = open(wigfilename, 'w')
        bw = pyBigWig.open(bwfilename,'w')

        bw.addHeader(list(sampleinfor.fregion.chrs_length.items()))

        kernel = smooth_kernel(kernellength)

        kernel_score = list()

        for w in sorted(kernel):

            kernel_score.append(kernel[w])

        #print('track type=wiggle_0 name="',sampleinfor.samplename,'" description="', sampleinfor.samplename,'"',sep='',file=wigio)

        for chromosome in sampleinfor.fregion.count_chr:

            chromosome = str(chromosome)

            #print('variableStep	chrom=', chromosome, sep='', file=wigio)

            # chrregion = chromosome+":"+str(1)+"-"+str(sampleinfor.fregion.chrs_length[chromosome])
            starts = []
            values = []

            for scare in range(0, int(sampleinfor.fregion.chrs_length[chromosome]/1000000)+1):

                startsite = scare * 1000000 + 1

                endsite = (scare + 1) * 1000000

                if endsite > sampleinfor.fregion.chrs_length[chromosome]:

                    endsite = sampleinfor.fregion.chrs_length[chromosome]

                # print (chromosome,startsite,endsite)

                # regionnow = chromosome+":" + str(startsite) + "-" + str(endsite)

                smoothedscore = regionsmooth(bamfile=bamfile, regionchromosome=str(chromosome),
                                             regionstart=startsite, regionend=endsite,
                                             chr_length=sampleinfor.fregion.chrs_length[chromosome],
                                             kernelsize=kernellength)

                # print (smoothedscore['score'])

                for site in sorted(smoothedscore['score'].keys()):

                    # print (site, smoothedscore['score'][site])

                    score = smoothedscore['score'][site]/samplenormalizedratio

                    if score > maxscore:

                        score = maxscore

                    score = round(score, 3)

                    starts.append(site)
                    values.append(score)

            bw.addEntries(chromosome, starts=starts, values=values,
                      span=1, step=1)

        bw.close()

                    #print (str(site)+"\t"+str(score), file=wigio)

        #wigio.close()

    except KeyboardInterrupt:

        print ("You cancelled the program!")

        sys.exit(1)