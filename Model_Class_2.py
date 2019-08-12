#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 10:54:47 2017

@author: thomas
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 17:51:38 2017

@author: thomas
"""

from numpy import (newaxis,array,ones,arange,stack,sqrt,cosh,
sinh,argwhere,sign,roll,isnan,linspace)
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt
from sympy import  N,symbols,nsolve


#import time


class Model_Layer_Media_array:
    
    def __init__(self, Model_Parameters):
        self.Model_Parameters = Model_Parameters
        self.step = 7
        self.Num_of_Layers = len(Model_Parameters)-1
        self.Half_Space_Parameters = Model_Parameters[len(Model_Parameters)-1]
        self.Layer_Parameters = Model_Parameters[0:len(Model_Parameters)-1]
        self.modes = []
        self.freqpermodes = []
        self.PandS = []
        self.stoneley = []
        self.v = []
        self.count = 0
        self.Det_Array_search_counter = 0
        self.Layer = arange(0,len(self.Layer_Parameters),1)
        
        """ Compute Alpha, Beta, and Stonely speeds """
        for x in self.Model_Parameters:
            alpha = sqrt(x[0]/x[4])
            beta = sqrt(x[3]/x[4])
            lPandS = [alpha,beta]
            self.PandS.append(lPandS)
            lstoneley = Model_Layer_Media_array.stone(x[0],x[1],x[2],x[3],x[4])
            self.stoneley.append(lstoneley)
            
    def stone(A,F,C,L,rho):
            """ computation of stoneley velocity """
            
            x = symbols('x')
            seq  = x + (((x-L)/(x-A))**.5)*((C*(x-A)+F**2)/(C*L)**.5)
            sol = nsolve(seq,0)
            csol = N(sol,n=12, chop=True)
            check = seq.subs(x,csol)
            velocity = (csol/rho)**.5
            return(csol,check,velocity)   
      
    def setstep(self,value):
        """set step size through velocity range """
        
        self.step = value
        
    def lp(self):
        """ print lists of layer properties """
        
        print('Layer Parametes (Includes Half Space)')
        print(self.Model_Parameters)
        print('PH and SH velocity (alpha/beta)')
        print(self.PandS)
        print('stoneley velocity')
        print(self.stoneley)
        
    def layer_value(self):
        """ Print values of layer model parameters to screen"""
        
        index = 1
        print('Layer n =   [c11,c13,c33,c44,rho,width]')
        print('#################################################################################')
        for x,y in zip(self.Model_Parameters,self.stoneley):
            print('Layer ',index,' = ',x,end = "")
            print(' #### alpha = ',sqrt(x[0]/x[4]),' #### beta = ',sqrt(x[3]/x[4]),end = "")
            print(' #### stoneley = ',y[2])
            index = index + 1
        print('Half Space = ',self.Half_Space_Parameters)
        print('Number of Layers = ',self.Num_of_Layers)
        return()
    
   
    
    
    def Ti_array(self, c, omega, LayerProp):
        
        #### Create wave number 2d array (omega X velocity) ####
        omega_array = array(omega)
        c_array = array(c)
        k_a = omega_array[:,newaxis]/c_array
        ########################################################
        
        #### Create arrays for each individual layer property ####
        c11,c13,c33,c44,rho,z = [], [], [], [], [], []
        for x in self.Layer:
            # Conditions of Layer i
            c11.append(LayerProp[x][0])# 1.98*10**10
            c13.append(LayerProp[x][1])
            c33.append(LayerProp[x][2])
            c44.append(LayerProp[x][3])  # .88*10**10
            rho.append(LayerProp[x][4]) # 2200
            z.append(LayerProp[x][5])  # 500
        
        c11_a = array(c11)
        c13_a = array(c13)
        c33_a = array(c33)
        c44_a = array(c44)
        rho_a = array(rho)
        z_a = array(z)
        ##########################################################
       
        #print('rho shape = ',rho_a.shape)
        omega_array_sqr = omega_array**2
        rho_mult_omega_sqr = omega_array_sqr[:,newaxis]*rho_a
        #print('rho_mult_omega_sqr',rho_mult_omega_sqr.shape)
        k_sqr = k_a**2
        I_k_a = ones(k_a.shape)
        #print('ika shape = ',I_k_a.shape)
        c11_3d = I_k_a * c11_a[:,newaxis,newaxis]
        c13_3d = I_k_a * c13_a[:,newaxis,newaxis]
        c33_3d = I_k_a * c33_a[:,newaxis,newaxis]
        c44_3d = I_k_a * c44_a[:,newaxis,newaxis]
        #print('param shape = ',c11_3d.shape)
        
        if c_array.ndim == 1: ## Check velocity array for initial dimension ##
            Io = ones(len(c_array)) ## Used to maintain proper shape ##
        else:
            Io = ones(len(c_array[0]))
        
        temp = []
        for x in rho_mult_omega_sqr:
            temp.append(Io*x[:,None])
        rmos_3d = stack((temp),axis=1)
        
        c33c44 = c33_3d*c44_3d
        c13plc44 = c13_3d+c44_3d
        c11ksqr = c11_3d*k_sqr
        
        M1 = c44_3d*(rmos_3d-c44_3d*k_sqr)+c33_3d*(rmos_3d-c11ksqr)+(c13plc44)**2*k_sqr
        M2 = (rmos_3d-c44_3d*k_sqr)*(rmos_3d-c11ksqr)
        
        radical = M1**2-4*c33c44*M2
        radical = radical.astype(complex)
        M1 = M1.astype(complex)
        c33c44 = c33c44.astype(complex)
        
        v1 = sqrt((-M1+sqrt(radical))/(2*c33c44))
        v3 = sqrt((-M1-sqrt(radical))/(2*c33c44))
        
        I_z = ones(c_array.shape)
        I_z = I_z.astype(complex)
        z_a2d = I_k_a * z_a[:,newaxis,newaxis] #z_a[:,newaxis]*I_z
        z_a2d = z_a2d.astype(complex)
        
        v1z2d = v1*z_a2d
        v3z2d = v3*z_a2d
        
        C1 = cosh(v1z2d)
        C3 = cosh(v3z2d)
        S1 = sinh(v1z2d)
        S3 = sinh(v3z2d)
        
        
        e1 = ((c13plc44)*k_a*v1)/(c11ksqr-rmos_3d-c44_3d*v1**2)
        e3 = ((c13plc44)*k_a*v3)/(c11ksqr-rmos_3d-c44_3d*v3**2)
        
        d1 = c33_3d*v1-c13_3d*k_a*e1
        d2 = c33_3d*v3-c13_3d*k_a*e3
        d3 = k_a+v1*e1
        d4 = k_a+v3*e3
        d5 = v1*e1-v3*e3
        d6 = v3*e1-v1*e3
        
        C1C3 = C1*C3
        S1S3 = S1*S3
        C3S1 = C3*S1
        C1S3 = C1*S3
        
        c33d6 = c33_3d*d6
        c44d5 = c44_3d*d5
        
        T16 = (1/(c33d6*c44d5))*(2*d1*d2*(C1C3-1)-(d1**2+d2**2)*S1S3)
        T11 = 1 - T16
        T12 = (1/(c33d6*d5))*(d1*d2*(d3+d4)*(C1C3-1)-(d2**2*d3+d1**2*d4)*S1S3)
        T13 = (1/(c44d5))*(-d2*C1S3+d1*C3S1)
        T14 = (1/(c33d6))*(d1*C1S3-d2*C3S1)
        T15 = (1/(c33d6*c44d5))*((e1*d2+e3*d1)*(1-C1C3)+(e1*d1+e3*d2)*S1S3)
    
        T21 = (1/(c33d6*d5)*(e1*e3*(d3+d4)*(1-C1C3)+(e1**2*d4+e3**2*d3)*S1S3))
        T22 = C1C3+T16
        T23 = (1/(c44d5))*(-e3*C1S3+e1*C3S1)
        T24 = (1/(c33d6))*(e1*C1S3-e3*C3S1)
        T25 = (1/(c33d6*c44d5))*(2*e1*e3*(1-C1C3)+(e1**2+e3**2)*S1S3)
        T26 = -T21
        T31 = T14
        T32 = (c44_3d/(c33d6))*(-d1*d4*C1S3+d2*d3*C3S1)
        T33 = C1C3
        T34 = ((c44d5)/(c33d6))*S1S3
        T35 = T24
        T36 = -T31
        T41 = T13
        T42 = (1/d5)*(d2*d3*C1S3-d1*d4*C3S1)
        T43 = ((c33d6)/(c44d5))*S1S3
        T44 = T33
        T45 = T23
        T46 = -T41
        T51 = T12
        T52 = (c44_3d/(c33d6*d5))*(2*d1*d2*d3*d4*(1-C1C3)+(d1**2*d4**2+d2**2*d3**2)*S1S3)
        T53 = T42
        T54 = T32
        T55 = T22
        T56 = -T12
        T61 = T16
        T62 = -T12
        T63 = -T13
        T64 = -T14
        T65 = -T15
        T66 = T11
        
        flag = 0
        p = [[] for x in range(6)]
        new = [[] for x in range(6)]
        
        
        for (a11,a12,a13,a14,a15,a16,
             a21,a22,a23,a24,a25,a26,
             a31,a32,a33,a34,a35,a36,
             a41,a42,a43,a44,a45,a46,
             a51,a52,a53,a54,a55,a56,
             a61,a62,a63,a64,a65,a66) \
        in zip(T11,T12,T13,T14,T15,T16,
               T21,T22,T23,T24,T25,T26,
               T31,T32,T33,T34,T35,T36,
               T41,T42,T43,T44,T45,T46,
               T51,T52,T53,T54,T55,T56,
               T61,T62,T63,T64,T65,T66):
            
            if flag == 0:
                #print(count)
                p[0] = a51 #a11
                p[1] = a52 #a12
                p[2] = a53 #a13
                p[3] = a54 #a14
                p[4] = a55 #a15
                p[5] = a56 #a16
                flag = 1
            else:
                #print(count)
                new[0] = p[0]*a11+p[1]*a21+p[2]*a31+p[3]*a41+p[4]*a51+p[5]*a61
                new[1] = p[0]*a12+p[1]*a22+p[2]*a32+p[3]*a42+p[4]*a52+p[5]*a62
                new[2] = p[0]*a13+p[1]*a23+p[2]*a33+p[3]*a43+p[4]*a53+p[5]*a63
                new[3] = p[0]*a14+p[1]*a24+p[2]*a34+p[3]*a44+p[4]*a54+p[5]*a64
                new[4] = p[0]*a15+p[1]*a25+p[2]*a35+p[3]*a45+p[4]*a55+p[5]*a65
                new[5] = p[0]*a16+p[1]*a26+p[2]*a36+p[3]*a46+p[4]*a56+p[5]*a66
                p[0] = new[0]
                p[1] = new[1]
                p[2] = new[2]
                p[3] = new[3]
                p[4] = new[4]
                p[5] = new[5]
                
        return(p)
    
    def Half_Array(self,c, omega):

    
        #k = omega/cL
        omega_array = array(omega)
        c_array = array(c)
        
        k_a = omega_array[:,newaxis]/c_array
        I_k_a = ones(k_a.shape)
        omega_2d = omega_array[:,newaxis]*I_k_a
        
        # Conditions of half space
        '''
        c11 = 10.985*10**10
        c44 = 4.16*10**10
        c13 = c11 - 2*c44
        c33 = c11
        rho = 2600
        '''
        c11 = self.Half_Space_Parameters[0]
        c13 = self.Half_Space_Parameters[1]
        c33 = self.Half_Space_Parameters[2]
        c44 = self.Half_Space_Parameters[3]
        rho = self.Half_Space_Parameters[4]
        
        
        ro2dsqr = rho*omega_2d**2
        ksqr = k_a**2
        
        M1 = c44*(ro2dsqr-c44*ksqr)+c33*(ro2dsqr-c11*ksqr)+(c13+c44)**2*ksqr
        M2 = (ro2dsqr-c44*ksqr)*(ro2dsqr-c11*ksqr)
        radicand = M1**2-4*c33*c44*M2
        
        v1 = sqrt((-M1+sqrt(radicand))/(2*c33*c44))
        v3 = sqrt((-M1-sqrt(radicand))/(2*c33*c44))
        
    
        # Calculated parameters from Conditions of Half Space
        e1 = ((c13+c44)*k_a*v1)/(c11*ksqr-ro2dsqr-c44*v1**2)
        # e2 = ((c13+c44)*k*v2)/(c11*k**2-rho*omega**2-c44*v2**2)
        e3 = ((c13+c44)*k_a*v3)/(c11*ksqr-ro2dsqr-c44*v3**2)
        # e4 = ((c13+c44)*k*v4)/(c11*k**2-rho*omega**2-c44*v4**2)
    
        d1 = c33*v1-c13*k_a*e1
        d2 = c33*v3-c13*k_a*e3
        d3 = k_a+v1*e1
        d4 = k_a+v3*e3
        d5 = v1*e1-v3*e3
        d6 = v3*e1-v1*e3
        
        # Generate Half Space vector
        Half = [[] for x in range(6)]
        Half[0] = d2-d1
        Half[1] = e3-e1
        Half[2] = -c44*d5
        Half[3] = -c33*d6
        Half[4] = c44*(d1*d4-d2*d3)
        Half[5] = d1-d2
           
        return(Half)
    
    def Det_Array(self,c, omega):
    
        p = self.Ti_array(c,omega, self.Layer_Parameters)
        h = self.Half_Array(c,omega)
        
        det = p[0]*h[0]+p[1]*h[1]+p[2]*h[2]+p[3]*h[3]+p[4]*h[4]+p[5]*h[5]
        
        return(det)
    
    def root_search(self,c, omega):
        """ initial root search through given velocity X omega space """
        det = self.Det_Array(c,omega)
        mask_det = ma.array(det, mask=isnan(det))
        
        zeroCross = []
        lpsigned = sign(mask_det)  # take sign of array values only
        # Shift array by 1 index subtract from unshifted array ...
        for x in lpsigned.real:
            signchange = ((roll(x, 1) - x) != 0).astype(int)
            signchange[0] = 0  # remove artifact from array shift last index to first
            # zeroCross = index of lp/c1 array where root exsists
            zeroCross.append(argwhere(signchange))
        omega_s =[]
        root_c = []
        mode_num = []
        for x,z in zip(zeroCross,omega):
            mode_index = 0
            for y in x:
                #print(mode_index)
                mode_num.append(mode_index)
                #print('?????',y)
                #root_c.append(arange(c[y[0]-1],c[y[0]]+1,1))
                root_c.append(linspace(c[y[0]-1],c[y[0]],10))
                omega_s.append(z)
                mode_index = mode_index +1
        df = pd.DataFrame(root_c)
        df2 = pd.DataFrame(mask_det)
        df.to_csv("file_0.csv")
        df2.to_csv("file_mdet_0.csv")
                
        return(det,zeroCross,lpsigned,root_c,omega_s,mode_num)
    
    def Det_Array_search(self, c_in, omega_in,  propertyList, precision):
        """ improving accuracy of found roots """
        self.Det_Array_search_counter += 1
        det = self.Det_Array(c_in,omega_in)
        mask_det = ma.array(det, mask=isnan(det))
        zeroCross = []
        lpsigned = sign(mask_det)  # take sign of array values only
        # Shift array by 1 index subtract from unshifted array ...
        #print('lpsigned = ',lpsigned)
        for x in lpsigned.real:
            signchange = ((roll(x, 1) - x) != 0).astype(int)
            signchange[0] = 0  # remove artifact from array shift last index to first
            # zeroCross = index of lp/c1 array where root exists
            zeroCross.append(argwhere(signchange))
        omega_s =[]
        root_c = []
        #print('zeroCross = ',zeroCross)
        #print('entering loop')
        roots = []
        for x,z,v in zip(zeroCross,omega_in,c_in):
            #print('x = ',x,'outer loop')
            for y in x:
                #print('###',v[y[0]],v[y[0]-1],v[y[0]]+precision, v[y[0]]-precision)
                roots.append(v[y[0]])
                #print('appending', arange(v[y[0]-1],v[y[0]]+precision,precision))
                #root_c.append(arange(v[y[0]-1],v[y[0]]+precision,precision))
                #print(linspace(v[y[0]]-precision,v[y[0]]+precision,num = 10))
                root_c.append(linspace(v[y[0]-1],v[y[0]]+precision,num = 10))
                omega_s.append(z)
        #print('###############################')  
        #print(Det_Array_search.counter)
        "Creating data frames for debugging inspection"
        #df = pd.DataFrame(root_c)
        #df2 = pd.DataFrame(mask_det)
        #df.to_csv("file_{}.csv".format(self.Det_Array_search_counter))
        #df2.to_csv("file_mdet_{}.csv".format(self.Det_Array_search_counter))
        if self.Det_Array_search_counter < 4:
            #print('return c',root_c)
            det,zeroCross,lpsigned,root_c,omega_s,lpsigned,mask_det,roots = self.Det_Array_search(root_c,omega_s,propertyList,precision/10)
        
        
        return(det,zeroCross,lpsigned,root_c,omega_s,lpsigned,mask_det,roots)
   
    def Create_Disp(self,c,omega):
        self.Det_Array_search_counter = 0
        
        det,zero,lp,nr,os,mode_n = self.root_search(c,omega)
        newdet,newzero,newlp,newnr,newos,lps,mam,roots = self.Det_Array_search(nr,os,self.Layer_Parameters,.1)
        
        velocity = [[] for x in range(max(mode_n)+1)]
        freq = [[] for x in range(max(mode_n)+1)]
        for x,y,z in zip(mode_n,newos,roots):
            velocity[x].append(z)
            freq[x].append(y)
            #print(x,y,z)
        
        for x,y in zip(freq,velocity):
            plt.plot(x,y,'o',markersize=2)
        return(velocity,freq,det,newdet,nr,newnr)
    
    
    
    
    
    