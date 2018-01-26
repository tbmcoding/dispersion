#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 13:06:32 2017

@author: thomas
"""
import Model_Data
import Model_Class
import Model_Class_2
from numpy import linspace,arange

#Iso1 = Model_Class.Model_Layer_Media(Model_Data.Iso1)
#Ti1 = Model_Class.Model_Layer_Media(Model_Data.Ti1)
#Ti_Alt = Model_Class.Model_Layer_Media(Model_Data.Ti_Alt)
#Ti_Alt50 = Model_Class.Model_Layer_Media(Model_Data.Ti_Alt50)
#Ti_Andy = Model_Class.Model_Layer_Media(Model_Data.Ti_And)
#Ti_And10 = Model_Class.Model_Layer_Media(Model_Data.Ti_And10)
#Ti_And6 = Model_Class.Model_Layer_Media(Model_Data.Ti_And6)
#Ti_And_M11 = Model_Class.Model_Layer_Media(Model_Data.Ti_And_M11)
#Ti_Andy_L6 = Model_Class.Model_Layer_Media(Model_Data.Ti_Andy_L6)
#T6 = Ti_Andy_L6
#New_Andy = Model_Class.Model_Layer_Media(Model_Data.New_Andy)

Iso2 = Model_Class_2.Model_Layer_Media_array(Model_Data.Iso2)
Iso3 = Model_Class_2.Model_Layer_Media_array(Model_Data.Iso3)
Ti_Andy = Model_Class_2.Model_Layer_Media_array(Model_Data.Ti_And)
chen = Model_Class_2.Model_Layer_Media_array(Model_Data.chen)
#chen2 = Model_Class_2.Model_Layer_Media_array(Model_Data.chen2)
haskell = Model_Class_2.Model_Layer_Media_array(Model_Data.haskell)
Ti1 = Model_Class_2.Model_Layer_Media_array(Model_Data.Ti1)
Ti_Alt = Model_Class_2.Model_Layer_Media_array(Model_Data.Ti_Alt)
Ti_Alt50 = Model_Class_2.Model_Layer_Media_array(Model_Data.Ti_Alt50)
Ti_Alt_s = Model_Class_2.Model_Layer_Media_array(Model_Data.Ti_Alt_s)
Ti_Alt_s50 = Model_Class_2.Model_Layer_Media_array(Model_Data.Ti_Alt_s50)

