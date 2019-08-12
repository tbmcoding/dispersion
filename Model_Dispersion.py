#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 13:06:32 2017

@author: thomas

This file sets up instances of Model_Class_2 with data sets from Model_Data

Model Class was an older version of the implementation. Model_Class_2 is operational but
still in need of refinement. Once this file is run instnces will be ready to generate dispersion curves.

Use Case

first set up search range

linspace sets search range of frequency ie. linspace(start,stop,number of data points)
arange sets search range for velocity values ie. arange(low velocity, high velocity,step size)
  
omega = linspace(0,100,num=75)
c = arange(1500,3600,7)

Then we call Create_Disp method to plot our dispersion curves and get dispersion data

w,e,r,t,y,u = Iso2.Create_Disp(c,omega)

the first two values returned are the velocity and frequency values, the others are iterations of determinant
values mostly for debugging purposes.

Good search ranges are model specific. I will try adding more search examples in the near future. 
"""
import Model_Data
import Model_Class_2
from numpy import linspace,arange

Iso2 = Model_Class_2.Model_Layer_Media_array(Model_Data.Iso2)
Iso3 = Model_Class_2.Model_Layer_Media_array(Model_Data.Iso3)
Ti_Andy = Model_Class_2.Model_Layer_Media_array(Model_Data.Ti_And)
chen = Model_Class_2.Model_Layer_Media_array(Model_Data.chen)
haskell = Model_Class_2.Model_Layer_Media_array(Model_Data.haskell)
Ti1 = Model_Class_2.Model_Layer_Media_array(Model_Data.Ti1)
Ti_Alt = Model_Class_2.Model_Layer_Media_array(Model_Data.Ti_Alt)
Ti_Alt_s = Model_Class_2.Model_Layer_Media_array(Model_Data.Ti_Alt_s)
Ti_Alt_s50 = Model_Class_2.Model_Layer_Media_array(Model_Data.Ti_Alt_s50)

