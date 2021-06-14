# Phase-correction functions

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
from skimage.restoration import unwrap_phase
from .array_functions import grating


def square_mask(xsize, ysize, center):
  '''
    
    Generates a square mask of dimensions (20px x 20px) 

    Parameters
    ----------
    xsize : int
        Total x-size of the pattern.
    ysize : int
        Total y-size of the pattern.
    center : list, np.array
        Center of the square, such as [x0, y0].

    Returns
    -------
    np.ndarray
        2D mask with 1 or 0 entries.

    '''  
  y,x = np.meshgrid(np.arange(ysize), np.arange(xsize))
  x_center, y_center = center[0], center[1]
  mask_x = np.abs(x-x_center) <= 10
  mask_y = np.abs(y-y_center) <= 10
  mask_x, mask_y = mask_x.astype(np.int), mask_y.astype(np.int)
  return mask_x * mask_y


def vary_patches(xsize, ysize, b, ref_center, i0, ref_only=False, int_map=0, phase_mask=0):
  '''
    
    Generates two (20x20)-patches in 225 possible configurations within a
    (300x300) region.

    Parameters
    ----------
    xsize : int
        Full x-size of the pattern.
    ysize : int
        Full y-size of the pattern.
    b : float
        Grating 'k vector'..
    ref_center : list, numpy.array
        Coordinates of the reference patch center.
        [x_center, y_center].
    i0 : int
        One of the 225 possible positions for the second patch.
    ref_only : bool, optional
        Choose between True (generates only the ref. patch) or
        False (generates the ref. patch + a second patch). 
        The default is False.
    phase_mask : numpy.ndarray, optional
        A 2D array, with dimensions (xsize, ysize), with the phase values. 
        The default is 0.

    Returns
    -------
    list 
        [2D pattern: pattern containing the two patches,
         x axis distance: [distance between ref_x and x central position,
                           distance between second_x and x central position],
         y axis distance: [distance between ref_y and y central position,
                           distance between second_y and y central position],
         d: distance between the two patches]

    '''  
  x0, y0 = ref_center[0], ref_center[1]
  x = np.arange(-7*20,7*21,20)
  y = np.arange(-7*20,7*21,20)  
#  x = np.delete(x, np.where(x==0))
#  y = np.delete(y, np.where(y==0))
  patch_pos = np.zeros((15*15,2))
  s = 0
  for i in range(len(x)):
    for j in range(len(y)):
      patch_pos[s] = np.array([x[i] + x0,
                               y[j] + y0])
      s += 1
  i0 = i0%(7*7)
  grating_full = grating(xsize, ysize, b, [xsize//2, ysize//2], 
                         int_map=int_map, phase_mask=phase_mask)
  center_mask = square_mask(xsize, ysize, ref_center)
  patch_mask = square_mask(xsize, ysize, patch_pos[i0])
  if ref_only==True:
    mask_total = center_mask
  else: 
    mask_total = center_mask + patch_mask
    
  xaxis_distance = [int(xsize/2-x0), int(xsize/2-patch_pos[i0][0])]
  yaxis_distance = [int(ysize/2-y0), int(ysize/2-patch_pos[i0][1])]
  
  d = np.sqrt((x0-patch_pos[i0][0])**2 + (y0-patch_pos[i0][1])**2)
  
  return [grating_full * mask_total, xaxis_distance, yaxis_distance, d]
  
        
def gauss_2d(xy, x0, y0, sigma, amp):
  ''' 
    Generates a 2D gaussian distribution.  

    Parameters
    ----------
    xy : np.ndarray
        Array containing: [x-values array, y-values array].
    x0 : float
        Center position (in x) of the 2D Gaussian function.
    y0 : float
        Center position (in y) of the 2D Gaussian function.
    sigma : float
        Standard deviation.
    amp : float
        Amplitude.

    Returns
    -------
    gaussian2d : np.ndarray
        (len(xy[0]), len(xy[1]))-array containing a 2D gaussian distribution.

  '''  
    
  y,x = np.meshgrid(np.arange(len(xy[1])), 
                    np.arange(len(xy[0])))    
  gaussian2d = amp*np.exp(-(x-x0)**2/(2*sigma**2) 
                          -(y-y0)**2/(2*sigma**2))
  gaussian2d = gaussian2d.ravel()
  return gaussian2d  


def interf_2d(xy_interp, ref_center, 
              sigma, d, focus, lamb, theta, varphi):
  '''
    
    Generates a 2D interference pattern bounded by a gaussian envelope.

    Parameters
    ----------
    xy_interp : np.ndarray
        Array containing [x-values array, y-values array].
    ref_center : float
        Center of the gaussian envelope.
    sigma : float
        Standard deviation of the gaussian envelope.
    d : float
        Distance between the two patterns in the last lens before the CCD.
    focus : int, float
        Focal length of the last lens before the CCD (in m!)
    lamb : float
        Beam's wavelength (in m!).
    theta : float
        Angle between the line that connects the two patterns on the DMD
        and the x-direction of the CCD.
    varphi : float
        Constant phase-difference 
        (usually employed on the setup's phase-correction).

    Returns
    -------
    np.ndarray
        1D array containing the .ravel() of the 2D pattern.

  '''  
  y, x = np.meshgrid(xy_interp[1].astype(np.float), xy_interp[0].astype(np.float))
  x0, y0 = float(ref_center[0]), float(ref_center[1])
  f3 = focus
  gauss_pkg = np.exp(-(x-x0)**2/(2*sigma**2) - (y-y0)**2/(2*sigma**2))
  lamb_new = lamb/(2*np.sin(d/(2*f3))) / 5.3e-6
  k_new = 2*np.pi/lamb_new
  interf = (1 + np.cos(k_new*(np.cos(theta)*(x-x0) + np.sin(theta)*(y-y0))
                             + varphi))
  
  total_interf = gauss_pkg*interf
  
  return total_interf.ravel()


def map_treatment(phase_map):     
  '''
    
    Performs phase-unwrapping and local-weight average on the phase-map.

    Parameters
    ----------
    phase_map : np.ndarray
        1D array containing the .ravel() of the phase-map obtained 
        from the patch-varying phase-correction protocol.

    Returns
    -------
    np.ndarray
        .ravel() of the phase-map with unwrapped phase and
        locally-weighted averaged values.

    '''  
  lsize = int(np.sqrt(len(phase_map)))
    
  pmap = unwrap_phase(np.reshape(phase_map, (lsize,lsize)))
  
  for i in range(lsize):     
    pmap[:,i] = lowess(pmap[:,i], np.arange(lsize), return_sorted=False)  

  for i in range(lsize):
    pmap[i,:] = lowess(pmap[i,:], np.arange(lsize), return_sorted=False)    
    
  pmap -= pmap[(lsize-1)//2,(lsize-1)//2]  
    
  return pmap.ravel()











def interf_2d_test(xy_interp, ref_center, 
              sigma, theta, lamb_new, varphi):
  y, x = np.meshgrid(xy_interp[1].astype(np.float), xy_interp[0].astype(np.float))
  x0, y0 = float(ref_center[0]), float(ref_center[1])
  gauss_pkg = np.exp(-(x-x0)**2/(2*sigma**2) - (y-y0)**2/(2*sigma**2))
  k_new = 2*np.pi/lamb_new
  interf = (1 + np.cos(k_new*(np.cos(theta)*(x-x0) + np.sin(theta)*(y-y0))
                             + varphi))
  
  total_interf = gauss_pkg*interf
  
  return total_interf.ravel()    


def map_treatment_test(phase_map):     
  pmap = unwrap_phase(np.reshape(phase_map, (31,31)))
  
  for i in range(len(pmap[0])):     
    pmap[:,i] = lowess(pmap[:,i], np.arange(len(pmap[0])), return_sorted=False)  

  for i in range(len(pmap[0])):
    pmap[i,:] = lowess(pmap[i,:], np.arange(len(pmap[0])), return_sorted=False)    
    
  pmap -= pmap[15,15]  
    
  return pmap.ravel()  

def vary_patches_test(xsize, ysize, b, ref_center, i0, ref_only=False, phase_mask=0):
  x0, y0 = ref_center[0], ref_center[1]
  x = np.arange(-15*20,15*21,20)
  y = np.arange(-15*20,15*21,20)  
#  x = np.delete(x, np.where(x==0))
#  y = np.delete(y, np.where(y==0))
  patch_pos = np.zeros((31*31,2))
  s = 0
  for i in range(len(x)):
    for j in range(len(y)):
      patch_pos[s] = np.array([x[i] + x0,
                               y[j] + y0])
      s += 1
  i0 = i0%(31*31)
  grating_full = grating(xsize, ysize, [xsize//2,ysize//2], b, phase_mask, 0)
  center_mask = square_mask(xsize, ysize, ref_center)
  patch_mask = square_mask(xsize, ysize, patch_pos[i0])
  if ref_only==True:
    mask_total = center_mask
  else: 
    mask_total = center_mask + patch_mask
    
  xaxis_distance = [int(xsize/2 - x0), int(xsize/2 - patch_pos[i0][0])]
  yaxis_distance = [int(ysize/2 - y0), int(ysize/2 - patch_pos[i0][1])]
  
  d = np.sqrt((x0-patch_pos[i0][0])**2 + (y0-patch_pos[i0][1])**2)
  
  return [grating_full * mask_total, xaxis_distance, yaxis_distance, d]


def vary_patches_int(xsize, ysize, b, ref_center, phase_mask=0):
  '''
    
    Generates a (20px) x (20px) patch and vary it along the DMD.

    Parameters
    ----------
    xsize : int
        Full x-size of the pattern.
    ysize : int
        Full y-size of the pattern.
    b : float
        Grating 'k vector'..
    ref_center : list, numpy.array
        Coordinates of the reference patch center.
        [x_center, y_center].
   phase_mask : numpy.ndarray, optional
        A 2D array, with dimensions (xsize, ysize), with the phase values. 
        The default is 0.

    Returns
    -------
    np.ndarray 
        2D pattern: pattern containing the two patches,
         
    '''  

  grating_full = grating(xsize, ysize, b, [xsize//2,ysize//2], phase_mask, 0)
  center_mask = square_mask_int(xsize, ysize, ref_center)
  mask_total = center_mask
    
  return grating_full * mask_total

def square_mask_int(xsize, ysize, center):
  '''
    
    Generates a square mask of dimensions (20px x 20px) 

    Parameters
    ----------
    xsize : int
        Total x-size of the pattern.
    ysize : int
        Total y-size of the pattern.
    center : list, np.array
        Center of the square, such as [x0, y0].

    Returns
    -------
    np.ndarray
        2D mask with 1 or 0 entries.

    '''  
  y,x = np.meshgrid(np.arange(ysize), np.arange(xsize))
  x_center, y_center = center[0], center[1]
  mask_x = np.abs(x-x_center) <= 10
  mask_y = np.abs(y-y_center) <= 10
  mask_x, mask_y = mask_x.astype(np.int), mask_y.astype(np.int)
  return mask_x * mask_y