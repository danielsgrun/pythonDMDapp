# -*- coding: utf-8 -*-

#Testing script.

#Created on Mon Jun 14 12:28:49 2021

#@author: Daniel


from pythonDMDapp import maps_dmd, applications

#maps_dmd.extract_int_map() # Extract intensity map
 
#maps_dmd.extract_phase_zernike() # Extract phase map

applications.move_bars(10, 5, use_phase_map=False, save_images=False) # Moving two-bar system

#applications.bragg_beams(40, 20, 30, use_phase_map=False, save_images=False) # Bragg beams

#applications.impurity(50, 20, 30, use_phase_map=False, save_images=False) # Spinning impurity
