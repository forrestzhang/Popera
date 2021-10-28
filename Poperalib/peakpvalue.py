from .Peak import *
from .Peakpvalue import *
from .sta import *
from .Hotspot import *
from .countreads import *
import sys
from .region import *
from multiprocessing import Pool
from .bgcount import *
import numpy as np

class KeyboardInterruptError(Exception):

    pass


def peakcount_nocontrol(hotspots, pvalue, bamfile, initiallength, minlength, nthreads, fregion):

    bpc = get_bpc(bamfile, hotspots, fregion.filted_region,nthreads=nthreads)

    pars = list()

    for hotspot in hotspots:

        # print (hotspot.chromosome,hotspot.start,hotspot.end,hotspot.hotspotid)

        whether_filter = False

        for i in range(int(hotspot.start/100)-1,int(hotspot.end/100)+1):

            if hotspot.chromosome in fregion.filted_region:
                if i in fregion.filted_region[hotspot.chromosome]:

                    whether_filter = True


        if not whether_filter:

            par = dict()

            par['pvalue'] = pvalue

            par['hotspot'] = hotspot

            par['bamfile'] = bamfile

            par['initiallength'] = initiallength

            #par['mutype'] = mutype

            par['minlength'] = minlength

            par['chrlength'] = fregion.chrs_length[hotspot.chromosome]

            par['bpc'] = bpc

            pars.append(par)

    pool=Pool(nthreads)

    chrpeaks = list()

    try:
        hotspotpeaks = pool.map(hotspotpeak_nocontrol, pars)

        for nowpeaks in hotspotpeaks:

            numberofpeak = len(nowpeaks)

            if numberofpeak == 0:

                pass

                # print ("no peak")

            else:

                for peak in nowpeaks:

                    # print ("peak",peak.parent,peak.chromosome,peak.start,peak.end,peak.peakid)

                    chrpeaks.append(peak)

            # print ()
        pool.close()

        return chrpeaks


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



def hotspotpeak_nocontrol(par):

    try:
        pvalue = par['pvalue']

        hotspot = par['hotspot']

        bamfile = par['bamfile']

        initiallength = par['initiallength']

        #mutype = par['mutype']
        
        minlength = par['minlength']

        chrlength=par['chrlength']

        bpc = par['bpc']

        start = hotspot.start

        end = hotspot.end

        hotspotlength = end - start + 1

        chromosome = hotspot.chromosome

        parentid = hotspot.hotspotid

        # regionstart = start
        #
        # regionend = end

        regionstart = start - hotspotlength

        regionend = end + hotspotlength

        peaks = list()

        if regionstart < 1:

            regionstart = 1

        if regionend>chrlength:

            regionend = chrlength

        #hotspotregio = chromosome + ':' + str(regionstart) + '-' + str(regionend)

        readscount = dict()

        readscount = dhsinglereadscounter(bamfile = bamfile, regionchromosome = chromosome, 
                                          regionstart = regionstart, regionend = regionend)

        totalreads = 0

        for i in range(regionstart,regionend):
            if i in readscount:
                # uniqcount = uniqcount + 1
                totalreads = totalreads + readscount[i]

        flank_lambda = (totalreads+0.0) /(regionend - regionstart +1) * initiallength

        #lambda_mu = (totalreads+0.0) /(regionend - regionstart +1) * initiallength

        avg_mu = bpc * initiallength

        #if mutype == 'avg':

        #    lambda_mu = avg_mu

        #if mutype == 'flank':

        #    lambda_mu = flank_lambda

        #if mutype == 'min':

        #    lambda_mu = min(flank_lambda, avg_mu)

        #if mutype == 'max':

        #    lambda_mu = max(flank_lambda, avg_mu)
        lambda_mu = min(flank_lambda, avg_mu)
        peaks = list()
        
        #if end-start+1 < initiallength:
        if end-start+1 < 1:

            pass

        else:

            region_site_pvalue = dict()

            region_point = dict()

            for window_start in range(start,end-initiallength+1):
                window_end = window_start + initiallength
                readsinwindow = 0
                for j in range(window_start,window_end+1):
                    if j in readscount:
                        readsinwindow = readsinwindow + readscount[j]
                now_pvalue = poissonpvalue(x=readsinwindow, mu =lambda_mu)
                # print (window_start,window_end,readsinwindow,lambda_mu,now_pvalue)
                region_site_pvalue[window_start] = now_pvalue

            #filter pvalue

            uniqpoint = list()

            for now_window_start in region_site_pvalue:
                if region_site_pvalue[now_window_start] < pvalue:
                    # print (now_window_start)
                    for now_site in range(now_window_start,now_window_start+initiallength):
                        region_point[now_site] = 1

            uniqpoint = list(region_point.keys())

            if region_point:

                #peaksregion = continueregion(points=uniqpoint, minlength=initiallength)
                peaksregion = continueregion(points=uniqpoint, minlength=minlength)

                peaks = peakspvalue(peaksregion,readscount,totalreads,regionstart,regionend,bpc,pvalue,parentid,chromosome)

                if len(peaks) > 1:

                    merge_peaks = merge_region(peaks,initiallength)

                    new_peaks = peakspvalue(merge_peaks,readscount,totalreads,regionstart,regionend,bpc,pvalue,parentid,chromosome)

                    if len(new_peaks) > 0:

                        peaks = new_peaks

        return peaks

    except KeyboardInterrupt:

        raise KeyboardInterruptError()

def peakspvalue(peaksregion,readscount,totalreads,regionstart,regionend,bpc,pvalue,parentid,chromosome):

    peaks = list()

    initid = 1

    for now_region in peaksregion:

        peaks_reads = 0

        start_site = now_region['start_site']

        end_site = now_region['end_site']

        for i in range(start_site, end_site+1):

            if i in readscount:

                peaks_reads = peaks_reads + readscount[i]

        regionlength = end_site - start_site +1

        region_flank_lambda = (totalreads+0.0) / (regionend - regionstart +1) * regionlength

        avg_region_lambda = bpc * regionlength

        #if mutype == 'avg':

        #    region_lambda_mu = avg_region_lambda

        #if mutype == 'flank':

        #    region_lambda_mu = region_flank_lambda

        #if mutype == 'min':

        #    region_lambda_mu = min(region_flank_lambda, avg_region_lambda)

        #if mutype == 'max':

        #    region_lambda_mu = max(region_flank_lambda, avg_region_lambda)
        region_lambda_mu = min(region_flank_lambda, avg_region_lambda)
        minsite = 0

        peakscore = 0

        peaksite = 1

        for i in range(start_site,end_site+1):
            #if i in region_site_pvalue:
            #    totalpvalue = totalpvalue + region_site_pvalue[i]
            #else:
            #    print ("can't find p in", i)

            if i in readscount:

                if readscount[i] > peakscore:

                    peakscore = readscount[i]

                    peaksite = i


        now_region_pvalue = poissonpvalue(x=peaks_reads, mu = region_lambda_mu)

        if now_region_pvalue < pvalue:

            peakid = parentid + '.' + str(initid)

            initid = initid + 1
            outpvalue = max(float(now_region_pvalue),2e-16)
            nowpeak = Peakpvalue(chromosome=chromosome, start=start_site, end=end_site,
                           pvalue=-np.log10(outpvalue), peakpoint=peaksite,
                           parent=parentid, peakid=peakid)

            peaks.append(nowpeak)
    return peaks

def merge_region(regions,initiallength):
    out_regions = list()
    out_start_site = regions[0].start
    out_end_site = regions[0].end
    for i in range(len(regions)-1):
        now_region = regions[i]
        next_region = regions[i+1]
        if next_region.start - now_region.end <= initiallength:
            #now_start_site = now_region['start_site']
            out_end_site = next_region.end
            x = 1
        else:
            next_start_site = next_region.start
            next_end_site = next_region.end
            if i == len(regions)-2:
                out_regions.append({'start_site':out_start_site,'end_site':out_end_site})
                out_regions.append({'start_site':next_start_site,'end_site':next_end_site})
                if i > 0:
                    x = 1
                else:
                    x = 0
            out_start_site = next_start_site
            out_end_site = next_end_site            
    if x==1:
        out_regions.append({'start_site':out_start_site,'end_site':out_end_site})
    return out_regions


def main():
    pass

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)