from numpy import *


def kde(z, w, xv):

    return sum(exp(-0.5*((z-xv)/w)**2)/sqrt(2*pi*w**2))


def smooth_kernel(length):

    if length % 2 == 0:

        length = length + 1

    bandwidth = (length - 1)/6.0

    one_kernel = dict()
    
    for pos in linspace(-(length-1)/2, (length-1)/2, length):
    
        one_kernel[int(pos)] = kde(pos, bandwidth, 0)
    
    return one_kernel


def smooth_kernel_adj(length, minscore):

    if length % 2 == 0:

        length = length + 1

    bandwidth = (length - 1)/6.0

    one_kernel = dict()

    for pos in linspace(-(length-1)/2, (length-1)/2, length):

        one_kernel[int(pos)] = kde(pos, bandwidth, 0)/minscore

    return one_kernel
