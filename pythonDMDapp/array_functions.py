## 2D-array function generator for DMD ##

import numpy as np
from scipy.special import eval_genlaguerre
  
def grating_disk(xsize, ysize, b, radius, center, int_map=[1], phase_mask=0, phi=0):
  '''
    
    Generates a grating disk (i.e., a grating limited by a circular region).

    Parameters
    ----------
    xsize : int
        Total x-size of the array.
    ysize : int
        Total y-size of the array.
    b : float
        Grating 'k vector'.
    radius : int
        Disk radius.
    center : list, np.array
        List or numpy.array containing [x_center, y_center] positions.
    phase_mask : numpy.array
        2D array, with dimensions (xsize, ysize), with the phase values. 
        The default is 0.
    phi : float
        Global phase-value to be applied on the whole pattern area. 
        The default is 0.

    Returns
    -------
    grating_disk : numpy.ndarray
        Disk region with a grating inside and zero in every position outside.

  '''

  y,x = np.meshgrid(np.arange(ysize), np.arange(xsize))
  x0, y0 = center[0], center[1]
  grating = 1./2*(1+np.cos(2*np.pi*(x-y)*b - phase_mask + phi))  
  grating = np.real(grating)
  
  if np.sum(int_map) == 0:
    thr = 0.5 
    grating[grating>=1.0*thr] = 1
    grating[grating<1.0*thr] = 0  

  else:
    profile = grating / int_map
    profile = profile/np.max(profile)  
    thr = 0.5 * (np.cos(np.arcsin(profile)) + 1) 
    grating[grating>=0.95*thr] = 1
    grating[grating<0.95*thr] = 0  

  mask = (x-x0)**2 + (y-y0)**2 < radius**2  
  mask = mask.astype(np.int)
  grating_disk = grating*mask
  return grating_disk 

def grating(xsize, ysize, b, center, int_map=5.8, phase_mask=0, phi=0):
  '''
    
    Generates a grating in the whole area of the (xsize x ysize)-array.

    Parameters
    ----------
    xsize : int
        Total x-size of the array.
    ysize : int
        Total y-size of the array.
    b : float
        Grating 'k vector'.
    center : list, np.array
        List or numpy.array containing [x_center, y_center] positions.
    phase_mask : numpy.array, optional
        A 2D array, with dimensions (xsize, ysize), with the phase values. 
        The default is 0.
    phi : float
        Global phase-value to be applied on the whole pattern area. 
        The default is 0.

    Returns
    -------
    grating : numpy.ndarray
        Generates a 2D array with a constant grating structure along its whole area. 

    '''
  y,x = np.meshgrid(np.arange(ysize), np.arange(xsize))
  x0, y0 = center[0], center[1]
  grating = 1./2*(1+np.cos(2*np.pi*(x-y)*b - phase_mask + phi))  
  grating = np.real(grating)

  if np.sum(int_map) == 0:
    thr = 0.5 ;
    grating[grating>=1.0*thr] = 1
    grating[grating<1.0*thr] = 0  

  else:
    profile = grating / int_map
    profile = profile/np.max(profile)  
    thr = 0.5 * (np.cos(np.arcsin(profile)) + 1) 
    grating[grating>=0.95*thr] = 1
    grating[grating<0.95*thr] = 0

  return grating

def prepare_images_DMD(images_array):
  '''
    
    Prepares the patterns to be used as input on the DMD.

    Parameters
    ----------
    images_array : np.ndarray
        (n, xsize, ysize)-array, where 'n' is the number of different patterns
        patterns. IMPORTANT: 'n' MUST MATCH THE 'num_img' from dmd_functions'
        'seq_alloc'. 
        
    Returns
    -------
    np.ndarray
        A numpy.ndarray containing all the patterns' .ravel() in sequence.

    '''
  imgSeq = []
  images_array[images_array>=0.5] = 255
  images_array[images_array!=255] = 0  
  
  if (np.shape(images_array)[1] != 1920 or np.shape(images_array)[2] != 1080):
    print("Please, input n XxY images as a single array (n, X, Y).")
  else:
    for i in range(0,np.shape(images_array)[0]):
      image_DMD = images_array[i].T  
      imgSeq.append(image_DMD.ravel())
  
  return np.concatenate(imgSeq)


def lg(xsize, ysize, b, center, pl, w, phase_mask=0):
  '''
    
    Generation of (l,p)-Laguerre-Gauss 2D patterns with grating structure.

    Parameters
    ----------
    xsize : int
        Total x-size of the array.
    ysize : int
        Total y-size of the array.
    b : float
        Grating 'k vector'.
    center : list, np.array
        List or numpy.array containing [x_center, y_center] positions.
    pl : tuple
        Tuple containing (l,p) of the LG patterns.
    w : float
        Can be tuned to generate a bigger / smaller 2D pattern.
    phase_mask : numpy.array, optional
        A 2D array, with dimensions (xsize, ysize), with the phase values. 
        The default is 0.

    Returns
    -------
    lg : np.ndarray
        (l,p)-Laguerre-Gauss grated pattern, with dimensions (xsize, ysize).

  '''
  
  p, l = pl[0], pl[1]
  y,x = np.meshgrid(np.arange(ysize), np.arange(xsize)) 
  x0,y0 = center[0], center[1]  
  env = grating(xsize, ysize, b, center, phase_mask=phase_mask 
                                                  + l*np.arctan2(x-x0,y-y0))
  r = np.sqrt((x-x0)**2 + (y-y0)**2)
  lg = r**abs(l) * eval_genlaguerre(p,l, 4*r**2/w**2) * np.exp(-r**2/w**2)
  lg = lg/np.max(lg)   
  lg *= env
  lg[lg>=0.5] = 1
  lg[lg<0.5] = 0
  return lg  


def bar(xsize, ysize, b, center, dims, theta=np.pi/4, int_map=5/8, phase_mask=0):
  '''
    
    Generates a bar structure on a grating environment.

    Parameters
    ----------
    xsize : int
        Total x-size of the array.
    ysize : int
        Total y-size of the array.
    b : float
        Grating 'k vector'.
    center : list, np.array
        List or numpy.array containing [x_center, y_center] positions.
    dims : list, np.array
        List or np.array containing [width, height] of the bar structure.
    phase_mask : numpy.array, optional
        A 2D array, with dimensions (xsize, ysize), with the phase values. 
        The default is 0..

    Returns
    -------
    np.ndarray
        Bar structure generated on a grating pattern, with dimensions (xsize, ysize).

    '''  
  y,x = np.meshgrid(np.arange(ysize), np.arange(xsize))
  x0, y0 = center[0], center[1]
  env = grating(xsize, ysize, b, center, int_map = int_map, phase_mask=phase_mask)
  sx, sy = dims[0], dims[1]
  mask_x = np.abs(np.cos(theta)*(x-x0) + np.sin(theta)*(y-y0)) <= sx/2
  mask_y = np.abs(np.sin(theta)*(x-x0) - np.cos(theta)*(y-y0)) <= sy/2
  return env * mask_x * mask_y


     