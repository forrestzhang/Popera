def effectregion(chrlength, windowsize, bw):

    """
        count effect region
        ===================--
                        --=================--

    """
    scare = int(chrlength/windowsize)

    efregions = dict()

    for i in range(0, scare+1):
        efregions[i] = dict()
        if i == 0:

            efregions[i]['ctstart'] = 1
            efregions[i]['ctend'] = int(windowsize + 1.5 * bw)
            efregions[i]['efstart'] = 1
            efregions[i]['efend'] = int(windowsize)
        elif i == scare:

            efregions[i]['ctstart'] = int(i * windowsize - 1.5 * bw)
            efregions[i]['ctend'] = int(chrlength)
            efregions[i]['efstart'] = int(i * windowsize + 1)
            efregions[i]['efend'] = int(chrlength)

        else:

            efregions[i]['ctstart'] = int(i * windowsize - 1.5 * bw)
            efregions[i]['ctend'] = int((i + 1) * windowsize + 1.5 * bw)
            efregions[i]['efstart'] = int(i * windowsize + 1)
            efregions[i]['efend'] = int((i + 1) * windowsize)

    return efregions


def continueregion(points, minlength=2):

    points.sort()

    start_index = 0

    end_index = 0

    continue_region = list()

    for index_now in range(1, len(points)):

        pre_index = index_now - 1

        if points[pre_index] + 1 == points[index_now]:

            if index_now == len(points) -1:

                if points[index_now] - points[start_index] + 1>= minlength :
                    #print (points[start_index], points[index_now])
                    region_now = dict()
                    region_now['start_site'] = points[start_index]
                    region_now['end_site'] = points[index_now]
                    continue_region.append(region_now)

            else:

                end_index = index_now

        else:

            if points[end_index] - points[start_index] + 1 >= minlength :

                #print (points[start_index], points[end_index])
                region_now = dict()
                region_now['start_site'] = points[start_index]
                region_now['end_site'] = points[end_index]
                continue_region.append(region_now)

            start_index = index_now

            end_index = index_now

    return continue_region


def windowregion(chr_length, site, windowsize, chromsome):

    windowstart = site - int(windowsize/2)

    windowend = site + int(windowsize/2)

    if windowstart < 1:

        windowstart = 1

    if windowend > chr_length:

        windowend = chr_length

    windowregion = chromsome+":"+str(windowstart) +'-'+str(windowend)

    return windowregion