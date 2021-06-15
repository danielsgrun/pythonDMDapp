# -*- coding: utf-8 -*-

#Applications with the DMD.

#Contains functions for: 
#    -- Bragg Spectroscopy: bragg_beams 
#    -- Vortex Generation: bars, impurity

#Created on Mon Jun 14 10:45:23 2021

#@author: Daniel SG

import os
import numpy as np
import matplotlib.pyplot as plt 
import time
from .phase_correction_functions import * 
from .hardware import DMD, Camera
from .array_functions import *
from .zernike_functions import mask_from_zernike_coeff

b = 7e-2

def bragg_beams(lamb_interf, v, act_time, 
                c = 0.2,
                cx = 0,
                cy = 0.2,
                use_phase_map = True,
                use_int_map = False,
                save_images = False):    
  '''
    

    Parameters
    ----------
    lamb_interf : int, float
        Desired wavelength (in micron!) for the interf. pattern.
    v : int, float
        Desired velocity (in micron/s !) for the interf. pattern.
    act_time : int
        Acting time (in seconds!) for the DMD to be on.
    c : int, float; optional
        Fine-tuning parameter. The default is 0.2.
    cx : int, float, optional
        Fine-tuning parameter, x-direction. The default is 0.
    cy : int, float, optional
        Fine-tuning parameter, y-direction. The default is 0.     
    use_phase_map : bool; optional
        The default is True.
    use_int_map : bool; optional
        The default is False.
    save_images : Bool, optional
        The default is False.

  '''
    
  Device = DMD() # initializes DMD Class
  Device.start_DMD() # starts DMD with Vialux .dll

  num_img = 9

  #ccd = Camera() # initializes CCD Class
  #ccd.start_ccd() # starts CCD 
  #ccd.set_exposure(0.01) # sets exposure time (in ms)

  xsize, ysize = Device.DMD.nSizeX, Device.DMD.nSizeY

  r = 20 # grated disk "radius"

  if use_phase_map == True:
    phase_map_bckp = np.load('mask_zernikes.npy')
  else:
    phase_map_bckp = np.zeros((xsize, ysize))  

  if use_int_map == True:
    int_map = np.load('intensity_map_2804.npy')
    int_map = (int_map+1)/np.max(int_map+1)
  else:
    int_map = [0] 

  image_2d = np.zeros((num_img,xsize,ysize))

  f = 75e3 # last-lens' focus (in um!)
  lamb = 401e-3 # beam's wavelength (in um!)
  
  gamma = num_img * v / lamb_interf # required refreshing rate on the DMD
  tau = int(1 / gamma * 1e6) # in us!! 
  
  print(tau)

  d = lamb * f / lamb_interf # distance (in um!) between patches to produce the
                           # desired wavelength on the interf. pattern

  d = d / (10.8) # conversion between mirror size and pixel (in preparation for the DMD)


  # position of each beam on the DMD
  pos1 = [int(xsize/2 - d/(2*np.sqrt(2.))) - c*20, 
          int(ysize/2 + d/(2*np.sqrt(2.))) - c*20]
  pos2 = [int(xsize/2 + d/(2*np.sqrt(2.))) - c*20,
          int(ysize/2 - d/(2*np.sqrt(2.))) - c*20]

  
  mask_x = mask_from_zernike_coeff(shape=[xsize, ysize],
                                   radius=400,
                                   center=[int(xsize/2-0*c*20), 
                                           int(ysize/2-0*c*20)],
                                   vec=[0,- cx*20, - cx*20])
  
  mask_y = mask_from_zernike_coeff(shape=[xsize, ysize],
                                   radius=400,
                                   center=[int(xsize/2-0*c*20), 
                                           int(ysize/2-0*c*20)],
                                   vec=[0,+ cy*20, - cy*20])
  
  for i in range(0,num_img): # preparing num_img images in sequence
                             # varying the phase of one of them 
    phi = 2*np.pi/num_img * i      

    beam_1 = grating_disk(xsize, ysize, b, r, pos1,
                        int_map=int_map,
                        phase_mask = 1*phase_map_bckp + mask_x + mask_y,
                        phi=phi)

    beam_2 = grating_disk(xsize, ysize, b, r, pos2,
                        int_map=int_map,
                        phase_mask = 1*phase_map_bckp + mask_x + mask_y, 
                        phi=0)

    image_2d[i] = beam_1 + beam_2

  
  # Preparing images for the DMD  
  imgSeq = imgSeq = prepare_images_DMD(image_2d)

  
  # Allocating images into DMD and running it.
  Device.seq_alloc(imgSeq, pictureTime=tau, num_img=num_img, run_DMD=True)
#  time.sleep(1.0)
#  plt.imshow(ccd.get_image(zoom=True).T, origin='lower')


  if save_images == True:  
    ccd=Camera()
    ccd.start_ccd()
    ccd.set_exposure(0.2)
    
    figuras = np.zeros((5*num_img, 60, 60))
    
    time.sleep(1)
    
    y,x = np.meshgrid(np.arange(60), np.arange(60))  
    
    z = (x-29)**2 + (y-30)**2
    
    for i in range(5*num_img):
      figuras[i] = ccd.get_image()
      time.sleep(tau * 1e-6)  
      
    for i in range(5*num_img):
      plt.figure()
      plt.imshow(figuras[i, :,:].T, cmap='inferno')
      plt.contourf(z.T, np.arange(0,50,5), cmap='bone_r', alpha=0.5)
      plt.title('lambda_interf = {0:3.2f} micron'.format(lamb_interf))
      plt.savefig('BraggBeams_{0:03d}.png'.format(i))
      
    ccd.close_ccd() 
    
#    time.sleep(act_time)
    
    Device.stop_DMD()

    plt.close('all')

    create_video = ['ffmpeg -y -r 30 ',
                  '-f image2 -s 1920x1080 -i ',
                  'BraggBeams_%03d.png ',
                  '-vcodec libx264 ',
                  '-crf 25  -pix_fmt yuv420p ',
                  'BraggBeams_video.mp4 ']
 
    create_video = ' '.join(create_video) 
    os.system(create_video)    

  else:
    time.sleep(act_time)  
    Device.stop_DMD()



# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


def move_bars(v, act_time, 
              c=0.2,
              cx=0.2,
              cy=0,
              use_phase_map = True, 
              use_int_map = False,
              save_images = False):
  '''
    

    Parameters
    ----------
    v : int, float
        Velocity of the bars.
    act_time : int, float
        Time interval (in seconds) to keep the DMD on.
    c : int, float, optional
        Fine-tuning parameter. The default is 0.2.
    cx : int, float, optional
        Fine-tuning parameter, x-direction. The default is 0.5.
    cy : int, float, optional
        Fine-tuning parameter, y-direction. The default is 0.    
    use_phase_map : Bool, optional
        Use (or not) the phase-map. The default is True.
    use_int_map : Bool, optional
        Use (or not) the intensity-map. The default is False.
    save_images : Bool, optional
        Save a sequence of images. Requires an Ueye CCD! 
        The default is False.

    Returns
    -------
    Two-bar structure acting on the fourier plane.

  '''
    
  
  if use_int_map == False:
    int_map = np.array([0])  
  else:
    int_map = np.load('intensity_map_0605.npy')
    int_map = (int_map+1)/np.max(int_map+1)
    
  if use_phase_map == False:
    phase_map = 0
  else:
    phase_map = np.load('mask_zernikes.npy') 
    
  Device = DMD()
  Device.start_DMD()

  num_img = 20

  xsize, ysize = Device.DMD.nSizeX, Device.DMD.nSizeY  
  
  dimensions = [15, 700]
  
  c1 = [xsize//2+dimensions[0]//2.-c*20, 
        ysize//2+dimensions[0]//2.-c*20]
  c2 = [xsize//2-dimensions[0]//2.-c*20, 
        ysize//2-dimensions[0]//2.-c*20]
    
  image_2d = np.zeros((num_img, xsize, ysize))

  mask_1 = mask_from_zernike_coeff(shape=[xsize, ysize],
                                   radius=400,
                                   center=[int(xsize/2-0*c*20), 
                                           int(ysize/2-0*c*20)],
                                   vec=[0,70,70])
  
  mask_x = mask_from_zernike_coeff(shape=[xsize, ysize],
                                   radius=400,
                                   center=[int(xsize/2-0*c*20), 
                                           int(ysize/2-0*c*20)],
                                   vec=[0,-35. - cx*20, -35. - cx*20])

  for i in range(num_img):
      
    mask_y = mask_from_zernike_coeff(shape=[xsize,ysize],
                                          radius=400,
                                          center=[int(xsize/2-0*c*20),
                                                  int(ysize/2-0*c*20)],
                                          vec=[0,20-2*i + cy*20, 
                                               2*i-20 - cy*20])
    
    bar1 = bar(xsize, ysize, b, center=c1, dims=dimensions, 
               int_map = int_map,
               phase_mask=mask_1+mask_x+mask_y+phase_map)
    
    bar2 = bar(xsize, ysize, b, center=c2, dims=dimensions, 
               int_map = int_map,
               phase_mask=mask_x+mask_y+phase_map)
    
    image_2d[i] = bar1+bar2
      
  print("Images ready. Allocating to DMD.")  
    
#  tau = 1e6  
  tau = int(5.993 / v * 1e6)

  imgSeq = prepare_images_DMD(image_2d)
  Device.seq_alloc(imgSeq, pictureTime=tau, num_img=num_img)

  if save_images == True:  
    ccd=Camera()
    ccd.start_ccd()
    ccd.set_exposure(0.05)
    
    figuras = np.zeros((num_img, 60, 60))
    
    time.sleep(2)
    
    for i in range(num_img):
      figuras[i] = ccd.get_image()
      time.sleep((tau * 1e-6))
    
    y,x = np.meshgrid(np.arange(60), np.arange(60))  
    
    z = (x-29)**2 + (y-30)**2
    
    for i in range(num_img):
      plt.figure()
      plt.imshow(figuras[i].T, cmap='inferno')
      plt.contourf(z.T, np.arange(0,50,5), cmap='bone_r', alpha=0.5)
      #plt.title('$\Lambda$ = {0:3.2f} $\mu m'.format(lamb_interf))
      plt.savefig('Bars_{0:03d}.png'.format(i))
      
    plt.close('all')  
      
    ccd.close_ccd() 
    time.sleep(act_time)
    Device.stop_DMD()

    create_video = ['ffmpeg -y -r 30 ',
                  '-f image2 -s 1920x1080 -i ',
                  'Bars_%03d.png ',
                  '-vcodec libx264 ',
                  '-crf 25  -pix_fmt yuv420p ',
                  'Bars_video.mp4 ']
 
    create_video = ' '.join(create_video) 
    os.system(create_video) 

  else:

    time.sleep(act_time)  
  
    Device.stop_DMD()



# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


def impurity(diameter, omega, act_time, 
                cx = 0.2, cy = 0,
                use_phase_map = True,
                use_int_map = False,
                save_images = False):    
  '''
    

    Parameters
    ----------
    diameter : int, float
        Desired diameter (in micron!) for the impurity circulation.
    omega : int, float
        Desired angular velocity (in rad/s !) for the impurity.
    act_time : int
        Acting time (in seconds!) for the DMD to be on.
    cx : int, float; optional
        x fine-tuning parameter. The default is 0.
    cy : int, float; optional
        y fine-tuning parameter. The default is 0.    
    use_phase_map : bool; optional
        The default is True.
    use_int_map : bool; optional
        The default is False.
    save_images : Bool, optional
        The default is False.

  '''
    
  Device = DMD() # initializes DMD Class
  Device.start_DMD() # starts DMD with Vialux .dll

  num_img = 15
  dtheta = 2*np.pi / num_img

  #ccd = Camera() # initializes CCD Class
  #ccd.start_ccd() # starts CCD 
  #ccd.set_exposure(0.01) # sets exposure time (in ms)

  xsize, ysize = Device.DMD.nSizeX, Device.DMD.nSizeY

  b = 7e-2 # grating "wave vector"
  r = 300 # grated disk "radius"

  if use_phase_map == True:
    phase_map_bckp = np.load('mask_zernikes.npy')
  else:
    phase_map_bckp = np.zeros((xsize, ysize))  

  if use_int_map == True:
    int_map = np.load('intensity_map_2804.npy')
    int_map = (int_map+1)/np.max(int_map+1)
  else:
    int_map = [0] 

  image_2d = np.zeros((num_img,xsize,ysize))
 
  # 
  pos = [int(xsize/2), 
          int(ysize/2)]
  

  amp = (diameter + 3.309)/(5.987) # based on fitting
  
  
  dt = int(2 * np.pi / (num_img * omega) * 1e6)

  print(dt)

  param1, param2 = [], []

  for i in range(0,num_img):
    param1.append(amp*np.cos(2*np.pi/num_img*i))
    param2.append(amp*np.sin(2*np.pi/num_img*i))
    
  param1, param2 = np.array(param1), np.array(param2)  
  
  for i in range(0,num_img): # preparing num_img images in sequence
                             # varying the phase of one of them    

    mask1 = mask_from_zernike_coeff(shape=[xsize,ysize],
                                          radius=r,
                                          center=[int(xsize/2),int(ysize/2)],
                                          vec=[0, param1[i] - 20*cx, -20*cx])
    
    mask2 = mask_from_zernike_coeff(shape=[xsize,ysize],
                                          radius=r,
                                          center=[int(xsize/2),int(ysize/2)],
                                          vec=[0, 20*cy, param2[i] - 20*cy])

    imp = grating_disk(xsize, ysize, b, r, pos,
                        int_map=int_map,
                        phase_mask = 1*phase_map_bckp + mask1 + mask2)

    image_2d[i] = imp

  
  # Preparing images for the DMD  
  imgSeq = prepare_images_DMD(image_2d)

  
  # Allocating images into DMD and running it.
  Device.seq_alloc(imgSeq, pictureTime=dt, num_img=num_img, run_DMD=True)
#  time.sleep(1.0)
#  plt.imshow(ccd.get_image(zoom=True).T, origin='lower')


  if save_images == True:  
    ccd=Camera()
    ccd.start_ccd()
    ccd.set_exposure(0.2)
    
    figuras = np.zeros((5*num_img, 60, 60))
    
    time.sleep(1)
    
    y,x = np.meshgrid(np.arange(60), np.arange(60))  
    
    z = (x-29)**2 + (y-30)**2
    
    for i in range(5*num_img):
      figuras[i] = ccd.get_image()
      time.sleep(dt * 1e-6)  
      
    for i in range(5*num_img):
      plt.figure()
      plt.imshow(figuras[i, :,:].T, cmap='inferno')
      plt.contourf(z.T, np.arange(0,50,5), cmap='bone', alpha=0.5)
      plt.title('Diameter = {0:3.2f} micron'.format(diameter))
      plt.savefig('Impurity_{0:03d}.png'.format(i))
      
    ccd.close_ccd() 
    
#    time.sleep(act_time)
    
    Device.stop_DMD()

    plt.close('all')

    create_video = ['ffmpeg -y -r 30 ',
                  '-f image2 -s 1920x1080 -i ',
                  'Impurity_%03d.png ',
                  '-vcodec libx264 ',
                  '-crf 25  -pix_fmt yuv420p ',
                  'Impurity_video.mp4 ']
 
    create_video = ' '.join(create_video) 
    os.system(create_video)    

  else:
    time.sleep(act_time)  
    Device.stop_DMD()
