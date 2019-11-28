

import sys
import logging
import os
from scipy.stats.mstats import *
from optparse import OptionParser
import pysam
from Poperalib.FRegion import *
from Poperalib.Hotspot import *
from Poperalib.hotspotscount import *
from Poperalib.bgcount import *
# from Poperalib.peakcount import *
from Poperalib.poperaio import *
from Poperalib.Sampleinfor import *
import re
from multiprocessing import Pool


def main():

    opt = opt_check(get_optparser())

    nocontrol(opt)

def readbed(bedfile):

    hotspots = []

    idnumber = 0

    samplename = bedfile.split('.')[0]

    with open(bedfile) as bed:

        for i in bed:

            inf = i.rstrip().split('\t')

            chrom = inf[0]

            start = int(inf[1])

            end = int(inf[2])

            idnumber += 1

            hotspotid = chrom + '.' + samplename +str(idnumber)

            hotspot = Hotspot(start=start, end=end, chromosome=chrom,
                                  hotspotid=hotspotid)

            hotspots.append(hotspot)

    bed.close()

    return hotspots

def nocontrol(opt):

    countchr = opt.countchr

    bedfilesstr = opt.bedfiles

    bedfiles = bedfilesstr.split(",")

    bamfilesstr = opt.datafiles

    bamfiles = bamfilesstr.split(',')

    samplenamesstr = opt.samplenames

    samplenames = samplenamesstr.split(",")

    minlength = opt.minlength

    nthreads = opt.nthreads


    premertedsite = dict()

    sampleinfors = list()

    lengthofdatafile = len(bamfiles)

    for i in range(0, lengthofdatafile):

        datafile = bamfiles[i]

        bedfile = bedfiles[i]

        samplename = samplenames[i]

        fregion = FRegion(bamfile=datafile,  countchr=countchr, nthreads=nthreads)

        sampleinfor = Sampleinfor(samplename=samplename, datafile=datafile, fregion=fregion)

        sampleinfors.append(sampleinfor)

        hotspots = readbed(bedfile)

        singlehotsportcountwrite(sampleinfor,hotspots,nthreads)

        for hotspot in hotspots:

            if hotspot.chromosome in premertedsite:

                for nowsite in range(hotspot.start, hotspot.end+1):

                    if nowsite in premertedsite[hotspot.chromosome]:

                        pass

                    else:

                        premertedsite[hotspot.chromosome][nowsite] = 1
            else:

                premertedsite[hotspot.chromosome] = dict()

                for nowsite in range(hotspot.start, hotspot.end+1):

                    if nowsite in premertedsite[hotspot.chromosome]:

                        pass

                    else:

                        premertedsite[hotspot.chromosome][nowsite] = 1

    mergedhotspots = list()

    idnumber = 0

    for chromosome in premertedsite:

        hotspots_list = list(premertedsite[chromosome].keys())

        chrhostspots = continueregion(points=hotspots_list, minlength=minlength)

        for hotspotsnow in chrhostspots:

            start_site = hotspotsnow['start_site']

            end_site = hotspotsnow['end_site']

            idnumber = idnumber + 1

            hotspotid = chromosome + '.' + "merged" +str(idnumber)

            hotspot = Hotspot(start=start_site, end=end_site, chromosome=chromosome,
                                  hotspotid=hotspotid)

            mergedhotspots.append(hotspot)

    hotspotswriter(hotspots=mergedhotspots, samplename='merged_DHS_all')

    mergedhotsportscountwrite(sampleinfors=sampleinfors, mergedhotspots=mergedhotspots, nthreads=nthreads)

def get_optparser():

    usage = """usage: %prog <-d datafile> [-n name] [options]
    Example %prog -i dh_sample1.bed,dh_sample2.bed -d dh_sample1.bam,dh_sample2.bam -n sample1,sample2
    """

    description = "%prog DNase I hypersensitive site normalized reads count calculation"

    poperaopt = OptionParser(version="%prog 0.03", description=description, usage=usage, add_help_option=False)

    poperaopt.add_option("-h", "--help", action="help", help="show this help message and exit.")

    poperaopt.add_option("-d", "--data", dest="datafiles", type="string", help='data file, should be sorted bam format, example=DH_sample1.bam,DH_sample2.bam')

    poperaopt.add_option("-n", "--name", dest="samplenames", help="NH sample name default=DH_sample1,DH_sample2", type="string" , default="DH_sample")

    poperaopt.add_option("-i", "--bed", dest="bedfiles", type="string", help="bed file, example=DH_sample1.bed,DH_sample2.bed")

    poperaopt.add_option("-l", "--minlength", dest="minlength", type="int", help="minimum length of merged hot spots, default=5", default= 50)

    poperaopt.add_option("--threads", dest="nthreads", type="int", help="threads number or cpu number, default=4", default=4)


    return poperaopt


def opt_check(poperaopt):

    (opt, args) = poperaopt.parse_args()

    if not opt.datafiles:

        logging.error("you need input a bam file, '-d dh_sample1.bam,dh_sample2.bam'")

        poperaopt.print_help()

        sys.exit(1)

    datafiles = opt.datafiles.split(",")

    for datafile in datafiles:

        if not os.path.isfile(datafile):

            logging.error("No such file: %s" % opt.datafile)

            sys.exit(1)

        dataindexfile = datafile + '.bai'

        if not os.path.isfile(dataindexfile):

            logging.error("Missing bam index file: %s" % dataindexfile)

            sys.exit(1)

    samplenames = opt.samplenames.split(",")


    if not len(samplenames) == len(datafiles):

        logging.error("The number is different between samplenames and datafiles")

        sys.exit(1)


    if not (opt.nthreads > 0):

        logging.error("threads number should >=1")

        poperaopt.print_help()

        sys.exit(1)

    opt.countchr = list()

    samfile = pysam.Samfile(datafiles[0])

    sam_ref = samfile.references

    for i in sam_ref:

        opt.countchr.append(str(i))

    return opt

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)

