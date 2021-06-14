# ----------------------------------------------------------------------------
# Code completely taken from: https://www.wavefrontshaping.net/
# ----------------------------------------------------------------------------

import numpy as np
from aotools.functions import phaseFromZernikes
  
def get_disk_mask(shape, radius, center = None):
    '''
    Generate a binary mask with value 1 inside a disk, 0 elsewhere
    :param shape: list of integer, shape of the returned array
    :radius: integer, radius of the disk
    :center: list of integers, position of the center
    :return: numpy array, the resulting binary mask
    '''
    if not center:
        center = (shape[0]//2,shape[1]//2)
    X,Y = np.meshgrid(np.arange(shape[0]),np.arange(shape[1]))
    mask = (Y-center[0])**2+(X-center[1])**2 < radius**2
    return mask.astype(np.int)
  
def mask_from_zernike_coeff(shape, radius, center, vec):
    '''
    Generate a complex phase mask from a vector containting the coefficient of the first Zernike polynoms.
    :param DMD_resolution: list of integers, contains the resolution of the DMD, e.g. [1920,1200]
    :param: integer, radius of the illumination disk on the DMD
    :center: list of integers, contains the position of the center of the illumination disk
    :center: list of float, the coefficient of the first Zernike polynoms
    '''
    # Generate a complex phase mask from the coefficients
    zern_mask = phaseFromZernikes(vec,2*radius)
    # We want the amplitude to be 0 outside the disk, we fist generate a binary disk mask
    amp_mask = get_disk_mask([2*radius]*2,radius)
    # put the Zernik mask at the right position and multiply by the disk mask
    mask = np.zeros(shape = shape)
    mask[center[0]-radius:center[0]+radius,
         center[1]-radius:center[1]+radius] = zern_mask*amp_mask
    return mask
