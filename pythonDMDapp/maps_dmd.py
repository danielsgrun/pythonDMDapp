# -*- coding: utf-8 -*-

#Created on Mon Jun 14 09:48:31 2021

#author: Daniel SG

#Zernike pol. code based on (and partially taken from): https://www.wavefrontshaping.net/


import numpy as np
import matplotlib.pyplot as plt 
import time
from tqdm import tqdm
from .array_functions import *
from .phase_correction_functions import *
from .hardware import DMD, Camera


def extract_phase_zernike(max_order=15, exposure=0.01, plot=True):
    '''
    Extract the DMD's phase-map.

    Parameters
    ----------
    max_order: int; optional
        Maximum order for the Polynomials. The default is 15. 
    exposure : int, float; optional
        Camera exposure time (in ms!). The default is 0.01.
    plot : Bool; optional
        Chooses to save the pictures each step. The default is True.

    Returns
    -------
    File "phase_map_bckp.npy" containing the adjusted phase-map.

    '''

    from .zernike_functions import mask_from_zernike_coeff

    num_img = 1  

    Device = DMD()
    Device.start_DMD()
    
    ccd = Camera()
    ccd.start_ccd()

    ccd.set_exposure(exposure)
  
    plot = True # choose between True and False

    xsize, ysize = Device.DMD.nSizeX, Device.DMD.nSizeY

    images_2d = np.zeros((num_img, xsize, ysize))

    ## Define the image!

    s = 0

    params = np.zeros((15))

    #int_map = np.load('intensity_map_0605.npy')
    int_map = [0]

    b = 7e-2

    gain = np.zeros((81))

    for i in range(3,15):
      print("Starting loop", i)
      for j in tqdm(range(81)):
    
        params[i] = -4.0 + 0.1*j
    #    params[i] = 8+0.2*j
    
        mask = mask_from_zernike_coeff(shape=[xsize,ysize],
                                          radius=400,
                                          center=[int(xsize/2),int(ysize/2)],
                                          vec=params)
        for k in range(num_img):
          images_2d[k] = grating_disk(xsize, ysize, b, radius=400, 
                                  center = [xsize//2, ysize//2],
                                  phase_mask=mask, int_map=[0])


        imgSeq = prepare_images_DMD(images_2d)

        Device.seq_alloc(imgSeq, num_img=num_img, pictureTime=1e5, run_DMD=True)
  
        time.sleep(0.5)
  
        ccd_frame = ccd.get_image(zoom=True)
        #cost_frame = cost(frame,
        #                  mask_radius = 8)
    
        np.save('frame_zernike', ccd_frame)
     
        frame = np.load('frame_zernike.npy')
    
        cost_frame = np.max(frame)
  
        plt.figure()
        plt.imshow(frame.T, cmap='inferno', origin='lower')
        #plt.title("Zernike coeff. 17 = {0:1.2f}".format(vals[j]))
        plt.colorbar()
    
  #    if j<p:
        plt.savefig("image_{0:04d}.png".format(s))
    
        plt.close('all')
  
        time.sleep(1)
  
        Device.free_seq()
  
        s += 1  
  
    #    print("Iteration:", s)
    #    print("Gain:", cost_frame)
        gain[j] = cost_frame
  
      p = np.where(gain==max(gain))[0][0]  
      params[i] = -4.0 + 0.1*p 
  
    #    imgs.append(frame)

    mask = mask_from_zernike_coeff(shape=[xsize,ysize],
                                          radius=400,
                                          center=[int(xsize/2),int(ysize/2)],
                                          vec=params)

    np.save('phase_map_bckp.npy', mask)

  
    Device.stop_DMD()  



# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------



def extract_int_map(exposure=0.1):
    
  '''
    Extracts the intensity-map of the beam arriving at the DMD.

    Parameters
    ----------
    exposure : int, float; optional
        Exposure time for the Camera (in ms!). The default is 0.1.

    Returns
    -------
    Creates "intensity_map_bckp.npy", containing the normalized intensity-map. 

  '''

  from scipy.interpolate import interp2d
  from scipy.optimize import curve_fit

  Device = DMD() # initializes DMD Class
  Device.start_DMD() # starts DMD with Vialux .dll

  ccd = Camera() # initializes CCD Class
  ccd.start_ccd() # starts CCD 
  ccd.set_exposure(exposure) # sets exposure time (in ms)

  xsize, ysize = Device.DMD.nSizeX, Device.DMD.nSizeY

  n_submaps_x = (xsize-20)//20 + 1
  n_submaps_y = (ysize-20)//20 + 1

  x_central = np.linspace(10, xsize-10, n_submaps_x)
  y_central = np.linspace(10, ysize-10, n_submaps_y)

  int_map_small = np.zeros((n_submaps_x, n_submaps_y))

  b = 7e-2

  image_2d = np.zeros((1,xsize,ysize))

  for i in tqdm(range(0,n_submaps_x)):
 # or i in tqdm(range(35-5,61+5)):
  #for i in [9,10,11,12,13]:
    for j in range(0,n_submaps_y):
    #for j in range(16-5,36+5):
    
      patch_0 = vary_patches_int(xsize, ysize, b, [x_central[i], y_central[j]], 0)  
    
      image_2d[0] = patch_0
    
      imgSeq = prepare_images_DMD(image_2d)      
      Device.seq_alloc(imgSeq, pictureTime=2e5)
      time.sleep(0.5)
      ccd_ref = ccd.get_image(zoom=True)
      time.sleep(0.5)
      np.save('refFrame', ccd_ref)

      refFrame = np.load('refFrame.npy')
  #    refFrame = refFrame - np.min(refFrame)
    
  #    noiseFrame = np.load('noiseFrame.npy')
    
  #    refFrame = refFrame - noiseFrame
    
      xy = np.array([np.arange(len(refFrame[:])),
                     np.arange(len(refFrame[:]))])
    
      try:
        ref_fit = curve_fit(gauss_2d, xy, refFrame.ravel(), p0=(45, 30, 30, 100), 
                            bounds=([10,10,1,0],[70,70,80,255]))
      except RuntimeError: 
        ref_fit = [[0,0,1,0]]  
      x0, y0 = ref_fit[0][0], ref_fit[0][1]
      sigma = abs(ref_fit[0][2])
      intensity = ref_fit[0][3]
      #intensity = np.sum(refFrame)
    
      int_map_small[i,j] = intensity
    
  
  int_map_small /= np.max(int_map_small)    
  
  x = np.arange(n_submaps_x)
  y = np.arange(n_submaps_y)

  xn = np.linspace(0, n_submaps_x, 1920)
  yn = np.linspace(0, n_submaps_y, 1080)

  int_interp = interp2d(y,x,int_map_small)

  int_map = int_interp(yn, xn) 

  np.save('intensity_map_bckp.npy', int_map)
