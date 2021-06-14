#Created on Sat Jun 12 12:32:32 2021

#author: Daniel SG

#Hardware-software interface for DMD & CCD

#CCD-related code based on (and partially taken from): GitHub::Mallekin/ueye_example.py

#DMD-related code based on (and partially taken from): GitHub::wavefrontshaping/ALP4lib


from ALP4 import *  # Necessary to connect w/ the DMD
from pyueye import ueye # Necessary to connect w/ the CCD
import numpy as np


class DMD:
  def __init__(self):
    '''
      Class DMD.

      Returns
      -------
      None.

      '''
    pass
  
  def start_DMD(self):
    
    '''
    Starts the DMD with the DLL provided by ALP4.    
    '''
      
    try:
      self.DMD.Halt()
    except:
      print("The DMD is not currently running.")
    else:
      print("DMD running has stopped.")
      
    try:
      self.DMD.FreeSeq()
    except:
      print("The DMD sequence is not currently allocated.")
    else:
      print("DMD sequence has stopped.")
      
    try:
      self.DMD.Free()
    except:
      print("The DMD is not currently idle.")
    else:
      print("DMD is free.")  
      
    self.DMD = ALP4(version='4.3')
    self.DMD.Initialize()
    
  def stop_DMD(self):
    '''
      Stops the DMD. 
      Requires 'start_DMD' for working again.
    '''
    
    try:
      self.DMD.Halt()
    except:
      print("The DMD is not currently running.")
    else:
      print("DMD running has stopped.")
      
    try:
      self.DMD.FreeSeq()
    except:
      print("The DMD sequence is not currently allocated.")
    else:
      print("DMD sequence has stopped.")
      
    try:
      self.DMD.Free()
    except:
      print("The DMD is not currently idle.")
    else:
      print("DMD is free.")
      
  def seq_alloc(self, imgSeq, num_img=1, pictureTime=1e6, run_DMD=True):
    '''
      
      Allocates a sequence of (SizeX x SizeY)-patterns to the DMD memory.
      Possibly runs the DMD (if chosen).

      Parameters
      ----------
      imgSeq : 1D numpy.ndarray
          A 1D list containing all the patterns concatenated (in order!)
               and in .ravel() configuration.     
      num_img : int
          Number of (SizeX x SizeY)-patterns to be allocated to the DMD.
          The default is 1.
      pictureTime : float
          Duration of each pattern on the DMD (in microseconds!). 
          The default is 1e6.
      run_DMD : bool
          Run the DMD imediately after allocating the patterns (True or False).
          The default is True.
          
      Returns
      -------
      None.

    '''
      
    try:  
      self.DMD.SeqAlloc(nbImg=num_img, bitDepth=1)
      self.DMD.SeqPut(imgData = imgSeq)
      self.DMD.SetTiming(pictureTime=int(pictureTime))  
    except: 
      print("Please, start the DMD first with 'start_DMD()'")
    if run_DMD==True:
      self.DMD.Run()
    else:
      pass

  def free_seq(self, verbose=True):
    '''
    
      Frees the DMD from the previous array sequence (if any).
    
      Parameters
      ----------
      verbose : bool
          Print outputs while freeing the DMD (True or False).
          The default is True.
          
      Returns
      -------
      None.

    '''
    try:
      self.DMD.Halt()
    except:
      if verbose==True:  
        print("The DMD is not currently running.")
      else: pass
    else:
      if verbose==True:  
        print("DMD running has stopped.")
      else: pass
    
    try:
      self.DMD.FreeSeq()
    except:
      if verbose==True:  
        print("The DMD sequence is not currently allocated.")
      else: pass
    else:
      if verbose==True:  
        print("DMD sequence has stopped.")
      else: pass


#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------


class Camera:
  def __init__(self):
    '''
      Class Camera

      Returns
      -------
      None.

      '''  
    pass
      
  def start_ccd(self):
    '''
      Starts the Ueye CCD.

      Returns
      -------
      Success / error.

      '''  

#Variables
    self.hCam = ueye.HIDS(0)             #0: first available camera;  1-254: The camera with the specified camera ID
    sInfo = ueye.SENSORINFO()
    cInfo = ueye.CAMINFO()
    self.pcImageMemory = ueye.c_mem_p()
    self.MemID = ueye.int()
    rectAOI = ueye.IS_RECT()
    self.pitch = ueye.INT()
    self.nBitsPerPixel = ueye.INT(24)    #24: bits per pixel for color mode; take 8 bits per pixel for monochrome
    channels = 3                    #3: channels for color mode(RGB); take 1 channel for monochrome
    m_nColorMode = ueye.INT()		# Y8/RGB16/RGB24/REG32
    self.bytes_per_pixel = int(self.nBitsPerPixel / 8)
    
    error = 0

# Starts the driver and establishes the connection to the camera
    nRet = ueye.is_InitCamera(self.hCam, None)
    if nRet != ueye.IS_SUCCESS:
      error += 1

# Reads out the data hard-coded in the non-volatile camera memory and writes it to the data structure that cInfo points to
    nRet = ueye.is_GetCameraInfo(self.hCam, cInfo)
    if nRet != ueye.IS_SUCCESS:
      error += 1

# You can query additional information about the sensor type used in the camera
    nRet = ueye.is_GetSensorInfo(self.hCam, sInfo)
    if nRet != ueye.IS_SUCCESS:
      error += 1

    nRet = ueye.is_ResetToDefault(self.hCam)
    if nRet != ueye.IS_SUCCESS:
      error += 1
    
    if error != 0:
      print("Error. Please, close any UEye-related apps and restart Camera.")  

# Set display mode to DIB
    nRet = ueye.is_SetDisplayMode(self.hCam, ueye.IS_SET_DM_DIB)

# Set the right color mode
    if int.from_bytes(sInfo.nColorMode.value, byteorder='big') == ueye.IS_COLORMODE_BAYER:
    # setup the color depth to the current windows setting
      ueye.is_GetColorDepth(self.hCam, self.nBitsPerPixel, m_nColorMode)
      self.bytes_per_pixel = int(self.nBitsPerPixel / 8)
      print("IS_COLORMODE_BAYER: ", )
      print("\tm_nColorMode: \t\t", m_nColorMode)
      print("\tself.nBitsPerPixel: \t\t", self.nBitsPerPixel)
      print("\tself.bytes_per_pixel: \t\t", self.bytes_per_pixel)
      print()

    elif int.from_bytes(sInfo.nColorMode.value, byteorder='big') == ueye.IS_COLORMODE_CBYCRY:
    # for color camera models use RGB32 mode
      m_nColorMode = ueye.IS_CM_BGRA8_PACKED
      self.nBitsPerPixel = ueye.INT(32)
      self.bytes_per_pixel = int(self.nBitsPerPixel / 8)
      print("IS_COLORMODE_CBYCRY: ", )
      print("\tm_nColorMode: \t\t", m_nColorMode)
      print("\tself.nBitsPerPixel: \t\t", self.nBitsPerPixel)
      print("\tself.bytes_per_pixel: \t\t", self.bytes_per_pixel)
      print()

    elif int.from_bytes(sInfo.nColorMode.value, byteorder='big') == ueye.IS_COLORMODE_MONOCHROME:
    # for color camera models use RGB32 mode
      m_nColorMode = ueye.IS_CM_MONO8
      self.nBitsPerPixel = ueye.INT(8)
      self.bytes_per_pixel = int(self.nBitsPerPixel / 8)
      print("IS_COLORMODE_MONOCHROME: ", )
      print("\tm_nColorMode: \t\t", m_nColorMode)
      print("\tself.nBitsPerPixel: \t\t", self.nBitsPerPixel)
      print("\tself.bytes_per_pixel: \t\t", self.bytes_per_pixel)
      print()

    else:
    # for monochrome camera models use Y8 mode
      m_nColorMode = ueye.IS_CM_MONO8
      self.nBitsPerPixel = ueye.INT(8)
      self.bytes_per_pixel = int(self.nBitsPerPixel / 8)
      print("else")

# Can be used to set the size and position of an "area of interest"(AOI) within an image
    nRet = ueye.is_AOI(self.hCam, ueye.IS_AOI_IMAGE_GET_AOI, rectAOI, ueye.sizeof(rectAOI))
    if nRet != ueye.IS_SUCCESS:
      print("is_AOI ERROR")

    self.width = rectAOI.s32Width
    self.height = rectAOI.s32Height

# Prints out some information about the camera and the sensor
    print("Camera model:\t\t", sInfo.strSensorName.decode('utf-8'))
    print("Camera serial no.:\t", cInfo.SerNo.decode('utf-8'))
    print("Maximum image self.width:\t", self.width)
    print("Maximum image self.height:\t", self.height) 
    print()


# Allocates an image memory for an image having its dimensions defined by self.width and self.height and its color depth defined by self.nBitsPerPixel
    nRet = ueye.is_AllocImageMem(self.hCam, self.width, self.height, self.nBitsPerPixel, self.pcImageMemory, self.MemID)
    if nRet != ueye.IS_SUCCESS:
      print("is_AllocImageMem ERROR")
    else:
    # Makes the specified image memory the active memory
      nRet = ueye.is_SetImageMem(self.hCam, self.pcImageMemory, self.MemID)
      if nRet != ueye.IS_SUCCESS:
        print("is_SetImageMem ERROR")
      else:
        # Set the desired color mode
        nRet = ueye.is_SetColorMode(self.hCam, m_nColorMode)



# Activates the camera's live video mode (free run mode)
    nRet = ueye.is_CaptureVideo(self.hCam, ueye.IS_DONT_WAIT)
    if nRet != ueye.IS_SUCCESS:
      print("is_CaptureVideo ERROR")

# Enables the queue mode for existing image memory sequences
    nRet = ueye.is_InquireImageMem(self.hCam, self.pcImageMemory, self.MemID, self.width, self.height, self.nBitsPerPixel, self.pitch)
    if nRet != ueye.IS_SUCCESS:
      print("is_InquireImageMem ERROR")


# Continuous image display
#while(nRet == ueye.IS_SUCCESS):
  def get_image(self,zoom=True):
    '''
      Obtains an image from the CCD.  
    
      Parameters
      ----------
      zoom : bool, optional
          Zoom into a given rectangular boundary.

      Returns
      -------
      np.ndarray
          Resulting image in a 2D array.

    '''
    # In order to display the image in an OpenCV window we need to...
    # ...extract the data of our image memory
    array = ueye.get_data(self.pcImageMemory, 
                          self.width, 
                          self.height, 
                          self.nBitsPerPixel, 
                          self.pitch, 
                          copy=False)

    # self.bytes_per_pixel = int(self.nBitsPerPixel / 8)

    # ...reshape it in an numpy array...
    frame = np.reshape(array,(self.height.value, 
                            self.width.value, self.bytes_per_pixel))

    # ...resize the image by a half
    frame = frame = frame[:,:,0]
     
    if zoom==True:
    
     #frame = frame[int(len(frame[:,0])/2)-147:int(len(frame[:,0])/2)-47, 
     #               int(len(frame[0,:])/2)-70:int(len(frame[0,:])/2)+30]
     
     frame = frame[475:535, 
                   555:615]
     
    else:
      frame = frame
        
    return frame.T    

    
  def set_exposure(self, exposure):
    """
    Set the exposure.

    Returns
    =======
    exposure: number
    Real exposure, can be slightly different than the asked one.
    """
    new_exposure = ueye.c_double(exposure)
    ueye.is_Exposure(self.hCam,
                     ueye.IS_EXPOSURE_CMD_SET_EXPOSURE,
                     new_exposure, 8)
    
    return new_exposure

  def get_exposure(self):
    """
    Get the current exposure.

    Returns
    =======
    exposure: number
        Current exposure.
    """
    exposure = ueye.c_double()
    ueye.is_Exposure(self.hCam, ueye.IS_EXPOSURE_CMD_GET_EXPOSURE,
                     exposure,  8)
    return exposure  


  def close_ccd(self):
    '''
      De-activates the Ueye CCD.

      Returns
      -------
      None.

      '''  

# Releases an image memory that was allocated using is_AllocImageMem() and removes it from the driver management
    ueye.is_FreeImageMem(self.hCam, self.pcImageMemory, self.MemID)

# Disables the hCam camera handle and releases the data structures and memory areas taken up by the uEye camera
    ueye.is_ExitCamera(self.hCam)

    print()
    print("Ueye CCD is off.")
