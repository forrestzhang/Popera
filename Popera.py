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
import re


def main():

    opt = opt_check(get_optparser())

    nocontrol(opt)


def nocontrol(opt):

    countchr = opt.countchr

    bw = opt.bw

    datafile = opt.datafile

    minlength = opt.minlength

    pvalue = opt.pvalue

    wig = opt.wig

    nthreads = opt.nthreads

    initiallength = opt.initial

    samplename = opt.samplename

    gff = opt.gff

    windowsize = int(1e5)

    threshold = opt.threshold

    fregion = FRegion(bamfile=datafile,  countchr=countchr, nthreads=nthreads)

    hotspots = hotspotscount_nocontrol(bamfile=datafile, threshold=threshold, kernellength=bw,
                                        windowsize=windowsize, nthreads=nthreads, minlength=minlength,
                                        samplename=samplename, fregion=fregion, countchr=countchr)

    hotspotswriter(hotspots=hotspots, samplename=samplename)


def get_optparser():

    usage = """usage: %prog <-i datafile> [-n name] [options]
    Example %prog -i dh_sample1.bam -n sample1
    """

    description = "%prog DNase I hypersensitive site identification"

    poperaopt = OptionParser(version="%prog 0.03 20140519", description=description, usage=usage, add_help_option=False)

    poperaopt.add_option("-h", "--help", action="help", help="show this help message and exit.")

    poperaopt.add_option("-d", "--data", dest="datafile", type="string", help='data file, should be sorted bam format')

    # poperaopt.add_option("-c", "--control", dest="controlfile", type="string", help='control(input) file, should be sorted bam format', default="no")

    poperaopt.add_option("-n", "--name", dest="samplename", help="NH sample name default=NH_sample", type="string" , default="DH_sample")

    poperaopt.add_option("-b", "--bandwidth", dest="bw", type="int", help="kernel smooth band width, should >1, default=50", default=200)

    poperaopt.add_option("-t", "--threshold", dest="threshold", type="float", help="Hot spots threshold, default=4.0", default=5.0)

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

    if not opt.datafile:

        logging.error("you need input a bam file, '-d dh_sample1.bam'")

        poperaopt.print_help()

        sys.exit(1)

    if not os.path.isfile (opt.datafile):

        logging.error("No such file: %s" % opt.datafile)

        sys.exit(1)

    dataindexfile = opt.datafile + '.bai'

    if not os.path.isfile (dataindexfile):

        logging.error("Missing bam index file: %s" % dataindexfile)

        sys.exit(1)

    # if not opt.controlfile == "no":
    #
    #     if not os.path.isfile (opt.controlfile):
    #
    #         logging.error("No such file: %s" % opt.controlfile)
    #
    #         sys.exit(1)
    #
    #     controlindexfile = opt.controlfile + '.bai'
    #
    #     if not os.path.isfile (controlindexfile):
    #
    #         logging.error("Missing bam index file: %s" % controlindexfile)
    #
    #         sys.exit(1)
    #
    # else:
    #
    #     opt.controlfile = "no"

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

    samfile = pysam.Samfile(opt.datafile)

    sam_ref = samfile.references

    for i in sam_ref:

        opt.countchr.append(str(i))

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
    # k = 0
    #
    # digstart = re.compile('\-')
    #
    # for m in opt.countchr:
    #
    #     if digstart.match(m):
    #
    #         del opt.countchr[k]
    #
    #         print("skip chr:", m)
    #
    #     k = k + 1
    #
    # else:
    #
    #     pass

    return opt

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)

