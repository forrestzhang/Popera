def cEM_zip(testdata):

    cdef float sumzip = sum(testdata)

    cdef int lengthoflist = len(testdata)

    cdef float phat = 0.5

    cdef float phatpre = -1.0

    cdef float lhatpre = -1.0

    cdef float base

    cdef int i = 0

    cdef int j = 0

    cdef int n

    cdef float c

    zhat = []

    zerolist = []

    for i  from 0<=i<lengthoflist:

        if testdata[i] == 0:

            zerolist.append(i)

    for j  from 0<=j<lengthoflist:

        zhat.append(0)

    lhat = sumzip/lengthoflist

#     print (lhat)

    while (i<1000000):

        i = i + 1
        base = ((phat + (1-phat) * 2.718281828459045 **(0-lhat)))
        for n in zerolist:

            if testdata[n] == 0:

                zhat[n] = phat/base

        sumzhat = sum(zhat)

        c = (lengthoflist-sumzhat)

        lhat = sumzip/c

        phat = sumzhat/lengthoflist

#         print (i, lhat, phat)

        if (abs(lhat-lhatpre)<0.000001 and abs(phat-phatpre)<0.000001):
        #if (abs(lhat-lhatpre)<0.000001):

            break

        else:

            phatpre = phat

            lhatpre = lhat

    return (lhat, phat)