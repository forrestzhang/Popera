

import pysam


def dhsinglereadscounter(bamfile, regionchromosome, regionstart, regionend):

    samfile = pysam.Samfile(bamfile)

    readscount = dict()

    for aligned_read in samfile.fetch(reference=str(regionchromosome), start=regionstart, end=regionend):

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


def dhsinglereadscounterstrand(bamfile, regionchromosome, regionstart, regionend):

    samfile = pysam.Samfile(bamfile)

    readscount = dict()

    readscount['-'] = dict()

    readscount['+'] = dict()

    for aligned_read in samfile.fetch(reference=str(regionchromosome), start=regionstart, end=regionend):

        if aligned_read.is_reverse:

            site = aligned_read.aend

        else:

            site = aligned_read.pos

        site = site + 1

        if aligned_read.is_reverse:

            if site in readscount['-']:

                readscount['-'][site] = readscount['-'][site] + 1

            else:

                readscount['-'][site] = 1
        else:

            if site in readscount['+']:

                readscount['+'][site] = readscount['+'][site] + 1

            else:

                readscount['+'][site] = 1

    return readscount


def dhsinglewindowscarecounter(bamfile, regionchromosome, regionstart, regionend, windowsize):

    samfile = pysam.Samfile(bamfile)

    window_size_count = dict()

    window_size_count['readscount'] = dict()

    window_size_count['uniquesite'] = 0

    # chromosome, sesite = region.split(':')
    #
    # startsite, endsite = sesite.split('-')
    #
    # startsite = int(startsite)
    #
    # endsite = int(endsite)

    uniqsites = dict()

    for alignend_read in samfile.fetch(reference=str(regionchromosome), start=regionstart, end=regionend):

        if alignend_read.is_reverse:

            site = alignend_read.aend

        else:

            site = alignend_read.pos

        if (regionstart<=site<=regionend):

            if site in uniqsites:
                pass
            else:
                uniqsites[site] = 1

        siteint=int(site/windowsize)

        if (regionstart<=siteint<=regionend):

            if siteint in window_size_count['readscount']:

                window_size_count['readscount'][siteint] = window_size_count['readscount'][siteint] + 1

            else:

                window_size_count['readscount'][siteint] = 1

    totaluniq = len(uniqsites)

    window_size_count['uniquesite'] = totaluniq

    return window_size_count


def dhsingleregioncounter(bamfile, regionchromosome, regionstart, regionend):

    samfile = pysam.Samfile(bamfile)

    window_count = 0

    # chromosome, sesite = region.split(':')
    #
    # startsite, endsite = sesite.split('-')
    #
    # startsite = int(startsite)
    #
    # endsite = int(endsite)

    for alignend_read in samfile.fetch(reference=str(regionchromosome), start=regionstart, end=regionend):

        if alignend_read.is_reverse:

            site = alignend_read.aend

        else:

            site = alignend_read.pos

        siteint=int(site)

        if (regionstart<=siteint<=regionend):

            window_count = window_count + 1

    return window_count


def windowmidsitecounter(bamfile, regionchromosome, regionstart, regionend, windowsize, chr_length):

    window_count = dict()

    # chromosome, sesite = region.split(':')
    #
    # startsite, endsite = sesite.split('-')
    #
    # startsite = int(int(startsite)-windowsize/2)
    #
    # endsite = int(int(endsite)+windowsize/2)

    if regionstart < 1:

        regionstart = 1

    if regionend > chr_length:

        regionend = chr_length

    # resizeregion = chromosome + ":" + str(startsite) + "-" + str(endsite)

    regioncount = dhsinglereadscounter(bamfile=bamfile, regionchromosome=str(regionchromosome),
                                       regionstart=regionstart, regionend=regionend)

    halfwindow = int(windowsize/2)

    midsite = regionstart

    while (midsite <= regionend):

        windowcount = 0

        for i in range(0-halfwindow,halfwindow):

            nowsite = midsite + i

            if nowsite in regioncount:

                nowcount = regioncount[nowsite]

                windowcount = windowcount + nowcount

        window_count[midsite] = windowcount

        midsite = midsite + 1

    return window_count


def dhsinglereadsnormailzed(bamfile, regionchromosome, regionstart, regionend, ultratio):

    samfile = pysam.Samfile(bamfile)

    readscount = dict()

    for aligned_read in samfile.fetch(reference=str(regionchromosome), start=regionstart, end=regionend):

        if aligned_read.is_reverse:

            site = aligned_read.aend

        else:

            site = aligned_read.pos

        if site in readscount:

            readscount[site] = readscount[site] + ultratio

        else:

            readscount[site] = ultratio

    return readscount


def dhsinglereadsexcount(bamfile, regionchromosome, regionstart, regionend, exsize=5):

    samfile = pysam.Samfile(bamfile)

    readscount = dict()

    for aligned_read in samfile.fetch(reference=str(regionchromosome), start=regionstart, end=regionend):

        if aligned_read.is_reverse:

            site = aligned_read.aend

        else:

            site = aligned_read.pos

        site = site + 1

        for nowsite in range(site-exsize, site+exsize):

            if nowsite in readscount:

                readscount[nowsite] = readscount[nowsite] + 1

            else:

                readscount[nowsite] = 1

    return readscount