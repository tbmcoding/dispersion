# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 15:37:37 2017

@author: thomas
"""

import DeltaBuild
from numpy import arange, argwhere, sign, roll
import matplotlib.pyplot as plt
from scipy import optimize


def Dispof(Conditions, frequency, rangeofc, step):
    # Create array of determinant values ################
    lp = []  # list to hold determinant values XY - ST
    c1 = arange(rangeofc[0], rangeofc[1], step,)  # range of velocity calcs
    # Main loop to calculate values of det
    for x in c1:
        lp.append(DeltaBuild.createDet(x, frequency,
                                     len(Conditions), Conditions))
    '''
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(c1,lp)
    plt.plot(c1, lp)
    fig.savefig('disfig.pdf')
    
    print(lp)
    print('')
    '''
    # Determine zero crossings #####################
    lpsigned = sign(lp)  # take sign of array values only
    # Shift array by 1 index subtract from unshifted array ...
    signchange = ((roll(lpsigned, 1) - lpsigned) != 0).astype(int)
    signchange[0] = 0  # remove artifact from array shift last index to first
    # zeroCross = index of lp/c1 array where root exsists
    zeroCross = argwhere(signchange)
    # Extract the roots ####################
    # function to take in single argument for optimize.bisect
    def fbisect(y):
        return(DeltaBuild.createDet(y, frequency, len(Conditions), Conditions))

    roots = []
    for x in range(len(zeroCross)):
        # Set bounds of bisection
        a = c1[zeroCross[x]-1]
        b = c1[zeroCross[x]+1]
        # extract roots using optimize.bisect from scipy package
        roots.append(optimize.bisect(fbisect, a, b))

    return(roots)
    