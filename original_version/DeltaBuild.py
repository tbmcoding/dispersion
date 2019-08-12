# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 14:39:37 2017

@author: thomas
"""

from sympy import Matrix, cosh, sinh, cos, sin


def createTi(c, omega, Layer, LayerProp):
    # wave number / velocity pair
    k = omega/c
    # Conditions of Layer i
    C_11 = LayerProp[Layer][0]  # 1.98*10**10
    C_44 = LayerProp[Layer][1]  # .88*10**10
    rho = LayerProp[Layer][2]  # 2200
    d = LayerProp[Layer][3]  # 500
    # Calculated parameters from Conditions of Layer i
    alpha = ((C_11)/rho)**.5
    beta = (C_44/rho)**.5
    # matrix element variable terms
    t = 2-(c**2/beta**2)
    gamma = beta**2/c**2
    mu = rho*beta**2

    if c == alpha or c == beta:  # avoidance of division by zero problem
        # print('Trip')
        c = c+1

    # velocity relations to determine proper eigenfunctions (from Buchan 1995)
    if c < alpha:
        r = (1-c**2/alpha**2)**(1/2)
        C_alpha = cosh(k*r*d)
        S_alpha = sinh(k*r*d)
    else:  # c Must be greater than alpha
        r = 1j*(c**2/alpha**2-1)**(1/2)
        r3 = (c**2/alpha**2-1)**(1/2)
        C_alpha = cos(k*r3*d)
        S_alpha = 1j*sin(k*r3*d)

    if c > beta:
        s = 1j*(c**2/beta**2-1)**(1/2)
        s3 = (c**2/beta**2-1)**(1/2)
        C_beta = cos(k*s3*d)
        S_beta = 1j*sin(k*s3*d)
    else:  # c Must be less than beta
        s = (1-c**2/beta**2)**(1/2)
        C_beta = cosh(k*s*d)
        S_beta = sinh(k*s*d)

    # Matrix elemetn variable terms depending on conditionals
    Q0 = ((1/(r*s)+r*s))*S_alpha*S_beta
    Q1 = ((t/(r*s))+2*r*s)*S_alpha*S_beta
    Q2 = (((t**2)/(r*s))+(4*r*s))*S_alpha*S_beta
    Q3 = ((t**3)/(r*s)+8*r*s)*S_alpha*S_beta
    Q4 = ((t**4)/(r*s)+16*r*s)*S_alpha*S_beta

    # Computation of Matrix Ti elements
    T11 = gamma**2*(-4*t+(t**2+4)*C_alpha*C_beta-Q2)
    T12 = (gamma**2/mu)*((2+t)*(1-C_alpha*C_beta)+Q1)
    T13 = (gamma/mu)*((1/s)*C_alpha*S_beta-r*S_alpha*C_beta)
    T14 = (gamma/mu)*(s*C_alpha*S_beta-(1/r)*S_alpha*C_beta)
    T15 = -T12
    T16 = (gamma/mu)**2*(2*(1-C_alpha*C_beta)+Q0)
    T21 = mu*gamma**2*(-2*t*(t+2)*(1-C_alpha*C_beta)-Q3)
    T22 = 1+(C_alpha*C_beta)-T11
    T23 = gamma*((t/s)*C_alpha*S_beta-2*r*S_alpha*C_beta)
    T24 = gamma*(2*s*C_alpha*S_beta-(t/r)*S_alpha*C_beta)
    T25 = 1-T22
    T26 = T12
    T31 = gamma*mu*(4*s*C_alpha*S_beta-(t**2/r)*S_alpha*C_beta)
    T32 = -T24
    T33 = C_alpha*C_beta
    T34 = -(s/r)*S_alpha*S_beta
    T35 = T24
    T36 = -T14
    T41 = gamma*mu*((t**2/s)*C_alpha*S_beta-4*r*S_alpha*C_beta)
    T42 = -T23
    T43 = -(r/s)*S_alpha*S_beta
    T44 = T33
    T45 = T23
    T46 = -T13
    T51 = -T21
    T52 = T25
    T53 = -T23
    T54 = -T24
    T55 = T22
    T56 = -T12
    T61 = gamma**2*mu**2*(8*t**2*(1-C_alpha*C_beta)+Q4)
    T62 = T21
    T63 = -T41
    T64 = -T31
    T65 = -T21
    T66 = T11
    # Assemble Layer Matrix Ti
    T = Matrix([[T11, T12, T13, T14, T15, T16], [T21, T22, T23, T24, T25, T26],
                [T31, T32, T33, T34, T35, T36], [T41, T42, T43, T44, T45, T46],
                [T51, T52, T53, T54, T55, T56], [T61, T62, T63, T64, T65, T66]])
    # Return the assembled Matrix Ti
    return(T)


def createHS(c, omega, numofLayers, propertyList):

    # Conditions of half space
    C_11d =  10.985*10**10
    C_44d =  4.16*10**10
    rho_d = 2600
    alpha_d = (C_11d/rho_d)**.5
    beta_d = (C_44d/rho_d)**.5
    # Calculated parameters from Conditions of Half Space

    if c == alpha_d or c == beta_d:
        c = c+1

    # velocity relations to determine proper eigenfunctions (from Buchan 1995)
    if c < alpha_d:
        r = (1-c**2/alpha_d**2)**(1/2)

    else:
        r = 1j*(c**2/alpha_d**2-1)**(1/2)

    if c > beta_d:
        s = 1j*(c**2/beta_d**2-1)**(1/2)

    else:
        s = (1-c**2/beta_d**2)**(1/2)

    # element terms
    mu = rho_d*beta_d**2
    t = 2-c**2/beta_d**2
    # Generate Half Space Matrix (from Buchan 1995)
    V = Matrix([[1-r*s],
                [mu*(t-2*r*s)],
                [mu*s*(2-t)],
                [-mu*r*(2-t)],
                [-mu*(t-2*r*s)],
                [mu**2*(4*r*s-t**2)]])

    return(V)


def createProp(c, omega, numofLayers, propertyList):
    # Loop through each layer creating Ti and Muliply Ti * Ti-1
    for x in range(numofLayers):
        # createTi args (velocity, frequency, Layer, Layer Param)
        Ti = createTi(c, omega, x, propertyList)
        if x == 0:
            # add free surface conditions at beginning of calculation
            U = Matrix([[0, 0, 0, 0, 0, 1]])
            prop = U * Ti
        else:
            prop = prop * Ti
    return(prop)


def createDet(c, omega, numofLayers, propertyList):
    U = Matrix([[0, 0, 0, 0, 0, 1]])
    UT = createProp(c, omega, numofLayers, propertyList)
    V = createHS(c, omega, numofLayers, propertyList)
    UTV = UT * V
   
    
    '''
    XY = UTV[0]*UTV[3]
    ST = UTV[1]*UTV[2]
    '''
    
    return(UTV[0])
