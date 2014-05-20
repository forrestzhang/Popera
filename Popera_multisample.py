from __future__ import division
from __future__ import print_function
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


def main():

    opt = opt_check(get_optparser())

    nocontrol(opt)


def nocontrol(opt):

    countchr = opt.countchr

    bw = opt.bw

    datafilesstr = opt.datafiles

    datafiles = datafilesstr.split(",")

    minlength = opt.minlength

    pvalue = opt.pvalue

    wig = opt.wig

    nthreads = opt.nthreads

    initiallength = opt.initial

    samplenamesstr = opt.samplenames

    samplenames = samplenamesstr.split(",")

    windowsize = int(1e5)

    threshold = opt.threshold

    premertedsite = dict()

    sampleinfors = list()

    lengthofdatafile = len(datafiles)

    for i in range(0, lengthofdatafile):

        datafile = datafiles[i]

        samplename = samplenames[i]

        fregion = FRegion(bamfile=datafile,  countchr=countchr, nthreads=nthreads)

        hotspots = hotspotscount_nocontrol(bamfile=datafile, threshold=threshold, kernellength=bw,
                                            windowsize=windowsize, nthreads=nthreads, minlength=minlength,
                                            samplename=samplename, fregion=fregion, countchr=countchr)



        for hotspot in hotspots:

            if hotspot.chromosome in premertedsite:

                for nowsite in range(hotspot.start, hotspot.end):

                    if nowsite in premertedsite[hotspot.chromosome]:

                        pass

                    else:

                        premertedsite[hotspot.chromosome][nowsite] = 1
            else:

                premertedsite[hotspot.chromosome] = dict()

                for nowsite in range(hotspot.start, hotspot.end):

                    if nowsite in premertedsite[hotspot.chromosome]:

                        pass

                    else:

                        premertedsite[hotspot.chromosome][nowsite] = 1

        sampleinfor = Sampleinfor(samplename=samplename, datafile=datafile, fregion=fregion)

        sampleinfors.append(sampleinfor)

        hotspotswriter(hotspots=hotspots, samplename=samplename)

    mergedhotspots = list()

    idnumber = 0

    for chromosome in premertedsite:

        hotspots_list = premertedsite[chromosome].keys()

        chrhostspots = continueregion(points=hotspots_list, minlength=minlength)

        for hotspotsnow in chrhostspots:

            start_site = hotspotsnow['start_site']

            end_site = hotspotsnow['end_site']

            idnumber = idnumber + 1

            hotspotid = chromosome + '.' + "merged" +str(idnumber)

            hotspot = Hotspot(start=start_site, end=end_site, chromosome=chromosome,
                                  hotspotid=hotspotid)

            mergedhotspots.append(hotspot)

    mergedfilename = samplenames[0]

    for i in range(1,lengthofdatafile):

        mergedfilename = mergedfilename+"_"+samplenames[i]

    hotspotswriter(hotspots=mergedhotspots, samplename=mergedfilename)

    mergedhotsportscountwrite(sampleinfors=sampleinfors, mergedhotspots=mergedhotspots, nthreads=nthreads)

    if wig:

        wigwritte(sampleinfors=sampleinfors, kernellength=bw, nthreads=nthreads)

def get_optparser():

    usage = """usage: %prog <-d datafile> [-n name] [options]
    Example %prog -d dh_sample1.bam,dh_sample2.bam -n sample1,sample2
    """

    description = "%prog DNase I hypersensitive site identification"

    poperaopt = OptionParser(version="%prog 0.03 20140519", description=description, usage=usage, add_help_option=False)

    poperaopt.add_option("-h", "--help", action="help", help="show this help message and exit.")

    poperaopt.add_option("-d", "--data", dest="datafiles", type="string", help='data file, should be sorted bam format, example=DH_sample1.bam,DH_sample2.bam')

    # poperaopt.add_option("-c", "--control", dest="controlfile", type="string", help='control(input) file, should be sorted bam format', default="no")

    poperaopt.add_option("-n", "--name", dest="samplenames", help="NH sample name default=DH_sample1,DH_sample2", type="string" , default="DH_sample")

    poperaopt.add_option("-b", "--bandwidth", dest="bw", type="int", help="kernel smooth band width, should >1, default=50", default=200)

    poperaopt.add_option("-t", "--threshold", dest="threshold", type="float", help="Hot spots threshold, default=4.0", default=4.0)

    poperaopt.add_option("-l", "--minlength", dest="minlength", type="int", help="minimum length of hot spots, default=5", default= 50)

    poperaopt.add_option("-p", "--pavlue", dest="pvalue", type="float", help="p-value cutoff for peak identification, default=0.05",
                        default=0.01)

    poperaopt.add_option("-i", "--initial", dest="initial", type="int", help="Peak's initial length, >1 and <minlength, default=5", default=5)

    poperaopt.add_option("--threads", dest="nthreads", type="int", help="threads number or cpu number, default=4", default=4)

    poperaopt.add_option("-w", "--wig", action="store_true", help="whether out put wiggle file, default=False", default=False)

    poperaopt.add_option("-f","--fdr",action="store_true",help="using FDR instead p-value", default=False)

    poperaopt.add_option("-x", "--excludechr", dest="excludechr", help="Don't count those DHs, example='-x ChrM,ChrC'")

    poperaopt.add_option("-g", "--gff", action="store_true", help="whether out put gff file, default=False", default=False)

    # poperaopt.add_option("-j","--jobtype",dest="jobtype",type="string",help="job type, such as nhpaired or nhsingle")

    # poperaopt.add_option("-m","--maxinsert",dest="maxinsert",type="int",help="when you use paired library, please set the maxinsert size",default=80)

    poperaopt.add_option("--pe", dest="pe", action="store_true", help="paired-end reads or single-end reads, default=False (single end)", default=False)

    return poperaopt


def opt_check(poperaopt):

    (opt, args) = poperaopt.parse_args()

    if not opt.datafiles:

        logging.error("you need input a bam file, '-d dh_sample1.bam,dh_sample2.bam'")

        poperaopt.print_help()

        sys.exit(1)

    datafiles = opt.datafiles.split(",")

    if len(datafiles) < 2:

        logging.error("Need at least 2 datafile: %s" % opt.datafiles)

        sys.exit(1)

    for datafile in datafiles:

        if not os.path.isfile(datafile):

            logging.error("No such file: %s" % opt.datafile)

            sys.exit(1)

        dataindexfile = datafile + '.bai'

        if not os.path.isfile(dataindexfile):

            logging.error("Missing bam index file: %s" % dataindexfile)

            sys.exit(1)

    samplenames = opt.samplenames.split(",")

    if len(datafiles) < 2:

        logging.error("Need at least 2 samplename: %s" % opt.datafiles)

        sys.exit(1)

    if not len(samplenames) == len(datafiles):

        logging.error("The number is different between samplenames and datafiles")

        sys.exit(1)

    if not (opt.pvalue > 0 and opt.pvalue < 1):

        logging.error("pvalue should be a float between 0 and 1")

        poperaopt.print_help()

        sys.exit(1)

    if not (opt.bw > 1):

        logging.error("band width should be a int greater than 1")

        poperaopt.print_help()

        sys.exit(1)

    if not (opt.initial >= 1 and opt.initial<=opt.minlength):

        logging.error("peak initial length should be >1 and <minlength")

        poperaopt.print_help()

        sys.exit(1)

    if not (opt.nthreads > 0):

        logging.error("threads number should >=1")

        poperaopt.print_help()

        sys.exit(1)

    opt.countchr = list()

    samfile = pysam.Samfile(datafiles[0])

    sam_ref = samfile.references

    for i in sam_ref:

        opt.countchr.append(i)

    if (opt.excludechr):

        excludchr = opt.excludechr.split(',')

        for chri in excludchr:

            if not chri in sam_ref:

                print (chri,'not in the %s file' % opt.datafile)

                print ("try to selcet exclude Chr from", end =" : ")

                print (sam_ref, sep=",")

                poperaopt.print_help()

                sys.exit(1)

            else:

                j = 0

                for n in opt.countchr:

                    if chri == n:

                        del opt.countchr[j]

                    j = j + 1

    ### skip digital start chromosome
    k = 0

    digstart = re.compile('^\d')

    for m in opt.countchr:

        if digstart.match(m):

            del opt.countchr[k]

            print("skip chr:", m)

        k = k + 1

    else:

        pass

    return opt

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)

