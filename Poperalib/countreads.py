from __future__ import division
from __future__ import print_function
import pysam


def dhsinglereadscounter(bamfile, region):

    samfile = pysam.Samfile(bamfile)

    readscount = dict()

    for aligned_read in samfile.fetch(region=region):

        if aligned_read.is_reverse:

            site = aligned_read.aend

        else:

            site = aligned_read.pos

        site = site + 1

        if site in readscount:

            readscount[site] = readscount[site] + 1

        else:

            readscount[site] = 1

    return readscount


def dhsinglewindowscarecounter(bamfile, region, windowsize):

    samfile = pysam.Samfile(bamfile)

    window_size_count = dict()

    window_size_count['readscount'] = dict()

    window_size_count['uniquesite'] = 0

    chromosome, sesite = region.split(':')

    startsite, endsite = sesite.split('-')

    startsite = int(startsite)

    endsite = int(endsite)

    uniqsites = dict()

    for alignend_read in samfile.fetch(region=region):

        if alignend_read.is_reverse:

            site = alignend_read.aend

        else:

            site = alignend_read.pos

        if (startsite<=site<=endsite):

            if site in uniqsites:
                pass
            else:
                uniqsites[site] = 1

        siteint=int(site/windowsize)

        if (startsite<=siteint<=endsite):

            if siteint in window_size_count['readscount']:

                window_size_count['readscount'][siteint] = window_size_count['readscount'][siteint] + 1

            else:

                window_size_count['readscount'][siteint] = 1

    totaluniq = len(uniqsites)

    window_size_count['uniquesite'] = totaluniq

    return window_size_count


def dhsingleregioncounter(bamfile, region):

    samfile = pysam.Samfile(bamfile)

    window_count = 0

    chromosome, sesite = region.split(':')

    startsite, endsite = sesite.split('-')

    startsite = int(startsite)

    endsite = int(endsite)

    for alignend_read in samfile.fetch(region=region):

        if alignend_read.is_reverse:

            site = alignend_read.aend

        else:

            site = alignend_read.pos

        siteint=int(site)

        if (startsite<=siteint<=endsite):

            window_count = window_count + 1

    return window_count




def windowmidsitecounter(bamfile, region, windowsize, chr_length):

    window_count = dict()

    chromosome, sesite = region.split(':')

    startsite, endsite = sesite.split('-')

    startsite = int(int(startsite)-windowsize/2)

    endsite = int(int(endsite)+windowsize/2)

    if startsite < 1:

        startsite = 1

    if endsite > chr_length:

        endsite = chr_length

    resizeregion = chromosome + ":" + str(startsite) + "-" + str(endsite)

    regioncount = dhsinglereadscounter(bamfile=bamfile, region=resizeregion)

    halfwindow = int(windowsize/2)

    midsite = startsite

    while (midsite <= endsite):

        windowcount = 0

        for i in range(0-halfwindow,halfwindow):

            nowsite = midsite + i

            if nowsite in regioncount:

                nowcount = regioncount[nowsite]

                windowcount = windowcount + nowcount

        window_count[midsite] = windowcount

        midsite = midsite + 1

    return window_count


def dhsinglereadsnormailzed(bamfile, region, ultratio):

    samfile = pysam.Samfile(bamfile)

    readscount = dict()

    for aligned_read in samfile.fetch(region = region):

        if aligned_read.is_reverse:

            site = aligned_read.aend

        else:

            site = aligned_read.pos

        if site in readscount:

            readscount[site] = readscount[site] + ultratio

        else:

            readscount[site] = ultratio

    return readscount


