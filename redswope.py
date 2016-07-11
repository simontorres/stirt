#!/usr/bin/python3.4
#Simon Torres 2016-05-03
#pipeline for SWOPE image reduction.

# no shutter correction

import sys
#this is here to avoid import all the modules innecesarily
if "--help" in sys.argv or "-h" in sys.argv:
  print("\nHelp Page:\n")
  print("Options:")
  print("\t--help\tor -h\t\t: Print this text")
  print("\t--debug\tor -v\t\t: Print more detailed information")
  print("\t--shutter or -sc\t: NI : Apply shutter correction")
  print("\t--skip-overscan\t\t: Don't re-do overscan correction nor trim")
  print("\t--skip-bias\t\t: Don't re-do bias correction")
  print("\t--skip-linearity\t: Don't re-do linearity correction")
  print("\t--skip-flats\t\t: Don't re-do flats")
  print("\t--data-path <opt>\t: Where is the raw data")
  print("\t--proc-path <opt>\t: Where to write processed data")
  print("\t--single-chip\t\t: Chip number, default is all.")
  print("\t--mosaic\t\t: NI : Build Data mosaic")
  print("\t--clean\t\t\t: NI : Remove files created")
  print("\n\tNI : Not Implemented (yet)")
  print("\n\n\tContact: simon@lco.cl\n")
  sys.exit(0)
#if len(sys.argv) < 3:
  ##sys.exit("Mode of Use: %s < params.txt > < dir-list.txt > < -v or --debug or --force-new-calib  >"%sys.argv[0])
  #sys.exit("Mode of Use: %s <data dir> < -v or --debug or --force-new-calib  >"%sys.argv[0])
#else:
  #dir_list = sys.argv[1]
import os
import glob
#import pyfits as fits
import numpy as np
import scipy.stats as stats
import re
import subprocess
import time
import warnings
import astropy.stats as asst
from astropy.io import fits
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
#cython stuff
import pyximport; pyximport.install()
import correct
#from itertools import groupby
#from operator import itemgetter

data_path	= "./"
proc_data_path	= "./"
search_pattern	= "ccd????c3.fits"
patt_trim_overscan  = 't_'
patt_bias_corrected = 'b_'
patt_shutter_correc = 's_'
patt_linearity_corr = 'l_'
patt_flat_corrected = 'f_'
bias_list	= []
flat_list	= []
filter_list	= []
sci_list	= []
flat_sci_list	= []
all_list	= []
quadrants	= [1,2,3,4]
linear_c2	= [0.033,0.012,0.010,0.014666]
linear_c3	= [0.0057,0.0017,0.0014,0.0027]
linear_alpha	= [1.1150,1.0128,1.00000,1.0696]

my_name = sys.argv[0].split("/")[-1]

def make_dir_list(target):
  sorted_list	= sorted(glob.glob(target+search_pattern))
  return sorted_list

def group(input):
  #input must be a 1-D array
  output = []
  for i in input:
    if i not in output:
      output.append(i)
  return output

#process related functions

def sort_list(ilist):
  global bias_list
  global flat_list
  global sci_list
  global flat_sci_list
  global all_list
  global filter_list
  
  list_all	= []
  filters	= []
  debug("Finding OBJECT type and Science FILTER")
  for i in range(len(ilist)):
    image	= ilist[i]
    data,header= get_data_header(image)
    objt	= header["OBJECT"]
    filt	= header["FILTER"]
    if "INSTRUME" not in header:
      print_spacers("ERROR!!!")
      debug("This tool only process images of the new CCD camera in SWOPE")
      sys.exit("I'm sorry, I Can't process this data")
    list_all.append([image,objt,filt])
    filters.append(filt)
    print_progress(i,len(ilist))
  filter_list	= group(filters)
  flat_list	= [[] for i in range(len(filter_list))]
  sci_list	= [[] for i in range(len(filter_list))]
  #print filter_list
  #print list_all
  debug("Classifying Images According to Observation Type or Filter")
  for e in range(len(list_all)):
    ii = list_all[e]
    image_root	= ii[0][2:-7]
    index_filter = filter_list.index(ii[2])
    if ii[1] == "bias":
      bias_list.append(image_root)
      all_list.append(image_root)
      #print ii[0][2:-7]
    elif ii[1] == "sflat":
      flat_list[index_filter].append(image_root)
      flat_sci_list.append(image_root)
      all_list.append(image_root)
    else:
      sci_list[index_filter].append(image_root)
      flat_sci_list.append(image_root)
      all_list.append(image_root)
    print_progress(e,len(list_all))
  debug("Filters detected\t= %s "%len(filter_list))
  debug("Bias images detected\t= %s "%len(bias_list))
  for i in range(len(filter_list)):
    debug("Flat Images detected in Filter %s\t= %s "%(filter_list[i],len(flat_list[i])))
    debug("Science Images detected in Filter %s\t= %s "%(filter_list[i],len(sci_list[i])))
  return 
      

#GET DATA
def get_data_header(file_name):
  try:
    hdu0	= fits.open(file_name)
    scidata	= hdu0[0].data
    header	= hdu0[0].header
    scidata = scidata.byteswap().newbyteorder().astype('float64')
    return scidata,header
  except IOError as err:
    debug("I/O error (get_data_header): %s"%err)
  except TypeError as err:
    debug("TypeError (get_data_header): %s : %s"%(file_name,err))
 
 
 
 
#OVERSCAN
#TRIM
def overscan_trim(image_list):
  length = len(image_list) * len(quadrants)
  for k in range(len(quadrants)):
    e = quadrants[k]
    for i in range(len(image_list)):
      image_root = image_list[i]
      current= i + k*len(image_list)
      print_progress(current,length)
      image = data_path + image_root + "c%s.fits"%e
      data,header = get_data_header(image)
      x,y = data.shape
      for j in range(x):
        data[j] = data[j]- np.median(data[j][2050:])
      new_data = data[:2056,:2048]
      proc_iname = proc_data_path+ patt_trim_overscan + re.sub(data_path,"",image)
      fits.writeto(proc_iname,new_data,header,clobber=True)
  return True



#BIAS
def make_bias(image_list):
  length = len(image_list) * len(quadrants)
  for k in range(len(quadrants)):
    e = quadrants[k]
    stack = []
    for i in range(len(image_list)):
      current= i + k*len(image_list)
      print_progress(current,length)
      image = image_list[i]
      full_image_name = proc_data_path+ patt_trim_overscan +image+"c%s.fits"%e
      #print full_image_name
      data,header = get_data_header(full_image_name)
      #if header["OBJECT"] == "bias":
      stack.append(data)
      #print_progress(e,len(image_list))
    header["OBJECT"] = "mbias"
    master = np.dstack(stack)
    master_bias = np.median(master,axis=2)
    master_bias_fname = proc_data_path + "master_bias_c%s.fits"%e
    fits.writeto(master_bias_fname,master_bias,header,clobber=True)
    debug("Master Bias is %s"%master_bias_fname)
  return master_bias



#BIASCORR
def bias_correct(image_list,master_bias):
  length = len(image_list) * len(quadrants)
  for k in range(len(quadrants)):
    e = quadrants[k]
    master_bias_name  = "master_bias_c%s.fits"%e
    master_bias, head = get_data_header(master_bias_name)
    for i in range(len(image_list)):
      current= i + k*len(image_list)
      print_progress(current,length)
      image_root = image_list[i]
      image	   = proc_data_path + patt_trim_overscan+image_root+"c%s.fits"%e
      data,header = get_data_header(image)
      bias_corrected_data = data-master_bias
      header["COMMENT"] = ("bias corrected image \ Added by %s"%my_name)
      bias_corrected_fname = proc_data_path + patt_bias_corrected +image_root+"c%s.fits"%e
      fits.writeto(bias_corrected_fname,bias_corrected_data,header,clobber=True)
  return

#FLAT
def flat(master):
  master = np.dstack(master)
  #master = asst.sigma_clip(master,axis=2,sig=3,cenfunc=np.median)
  cflat  = np.median(master,axis=2)
  nflat  = cflat/np.median(cflat)
  return nflat

#FLATS
def make_flats(image_list,filter):
  if len(image_list) > 1:
    for e in quadrants:
      master_flat = []
      for i in range(len(image_list)):
        length = len(image_list) * len(quadrants)
        current= i + (e-1)*len(image_list)
        print_progress(current,length)
        image = proc_data_path + patt_linearity_corr + image_list[i] + "c%s.fits"%e
        data,header  =get_data_header(image)
        master_flat.append(data)
      master = flat(master_flat)
      header["COMMENT"] = ("MasterFlat / Added by %s"%my_name)
      flat_name = proc_data_path + "master_flat_%s_c%s.fits"%(filter,e)
      fits.writeto(flat_name,master,header,clobber=True)
      debug("Master flat is %s"%flat_name)
    return
  else:
    debug("Not enough flat images")
    return

def flat_correct(image_list,filter):
  if len(image_list) > 0:
    print_spacers("Flat Correction")
    length = len(image_list) * len(quadrants)
    for k in range(len(quadrants)):
      e = quadrants[k]
      master_flat_name = proc_data_path + "master_flat_%s_c%s.fits"%(filter,e)
      #print master_flat_name
      master_flat, head= get_data_header(master_flat_name)
      for i in range(len(image_list)):
        image = proc_data_path + patt_linearity_corr + image_list[i] + "c%s.fits"%e
        data,header = get_data_header(image)
        flat_data = data / master_flat
        header["COMMENT"] = ("Flat Corrected \ Added by %s"%my_name)
        flat_image = proc_data_path + patt_flat_corrected + image_list[i] + "c%s.fits"%e
        fits.writeto(flat_image,flat_data,header,clobber=True)
        current= i + k*len(image_list)
        print_progress(current,length)
    return
  else:
    debug("Not enough science images")
    return 

def mosaic(image_list):
  if len(image_list) > 0:
    print_spacers("Mosaic")
    length = len(image_list)
    for i in range(len(image_list)):
      mosaic = []
      m_header = None
      for e in range(1,5):
        image= proc_data_path + patt_flat_corrected + image_list[i] + "c%s.fits"%e
        data,header = get_data_header(image)
        mosaic.append(data)
        #print(image)
      m_header = header
      #print(len(mosaic))
      chip_1 = mosaic[0][::-1,::-1]
      #print(type(chip_1))
      x1,y1 = chip_1.shape
      chip_2 = mosaic[1][::-1,:]
      x2,y2 = chip_2.shape
      chip_3 = mosaic[2]
      x3,y3 = chip_3.shape
      chip_4 = mosaic[3][:,::-1]
      x4,y4 = chip_4.shape
      chip_21 = np.concatenate((chip_2,chip_1),axis=1)
      #print(chip_21.shape)
      chip_34 = np.concatenate((chip_3,chip_4),axis=1)
      #print(chip_34.shape)
      final   = np.concatenate((chip_34,chip_21),axis=0)
      #print(final.shape)
      #print(type(final))
      output_name = proc_data_path + image_list[i] + ".fits"
      header.remove("DATASEC")
      fits.writeto(output_name,final,header,clobber=True)
      print_progress(i,length)
      #plt.imshow(final, cmap='gray',vmin=90,vmax = 300)
      #plt.colorbar()
      #plt.show()
      #sys.exit(0)
      #mosaic_image = np.empty([x3+x4,y3+y2])
      
      
      #print(image_list[i])
  return
    


#def linear(pix_val,c2,c3,alpha):
  #corr_pix = pix_val * ( 1.0 + c2*(pix_val/32000.)+c3*pow(pix_val/32000.,2))*alpha
  #return corr_pix

def linearity_correction(image_list):
  if "--shutter" in sys.argv:
    #self_seach_pattern = patt_shutter_correc
    self_seach_pattern = patt_bias_corrected # remove once shutter is implemented and uncomment previous
  else:
    self_seach_pattern = patt_bias_corrected
  length = len(image_list) * len(quadrants)
  for k in range(len(quadrants)):
    e = quadrants[k]
    for i in range(len(image_list)):
      current= i + k*len(image_list)
      print_progress(current,length)
      image	= proc_data_path + self_seach_pattern + image_list[i] + "c%s.fits"%e
      #print image
      data,header= get_data_header(image)
      header["COMMENT"] = ("linearity corrected image \ Added by %s"%my_name)
      c2	= linear_c2[e-1]
      c3	= linear_c3[e-1]
      alpha	= linear_alpha[e-1]
      lc_data = correct.linear(data,c2,c3,alpha)
      # = data
      lc_image = proc_data_path + patt_linearity_corr + image_list[i] + "c%s.fits"%e
      #print lc_image
      fits.writeto(lc_image,data,header,clobber=True)
      #print_progress(i,len(image_list))
  return


def shutter_corr(image_list):
  debug("Shutter correction not implemented Yet")

#aesthetic and user comunication functions

def print_progress(current,total):
  if current == total:
    sys.stdout.write("Progress {:.2%}\n".format(1.0*current/total))
  else:
    sys.stdout.write("Progress {:.2%}\r".format(1.0*current/total))
  sys.stdout.flush()
  return

def print_spacers(message):
  rows, columns = os.popen('stty size', 'r').read().split()
  if len(message)%2 == 1:
    message	= message+" "
  bar_length	= int(columns)
  bar 	= "="*bar_length
  blanks	= bar_length - 2
  blank_bar	= "="+" "*blanks +"="
  space_length	= int((blanks-len(message))/2)
  message_bar	= "=" + " " * space_length + message + " " * space_length + "="
  print(bar)
  #print(blank_bar)
  print(message_bar)
  #print(blank_bar)
  print(bar)
  
def remove_by_pattern(prefix):
  pattern 	= proc_data_path + patt_trim_overscan + "*"
  list_to_delete= glob.glob(pattern)
  for i in range(len(list_to_delete)):
    #remove instruction
    print_progress(i,len(list_to_delete))

def clean_dir():
  if '--clean' in sys.argv:
    patts  = [patt_trim_overscan,patt_bias_corrected,patt_linearity_corr,patt_flat_corrected]
    for pat in patts:
      remove_by_pattern(pat)
    debug("cleaning")

def debug(message):
  if '-v' in sys.argv or '--debug' in sys.argv:
    print(message)
  return

#def print_help():
  #print("\nHelp Page:\n")
  #print("Options:")
  #print("\t--help\tor -h\t\t: Print this text")
  #print("\t--debug\tor -v\t\t: Print more detailed information")
  #print("\t--shutter or -sc\t\t: NI : Apply shutter correction")
  #print("\t--skip-overscan\t\t: Don't do overscan correction nor trim")
  #print("\t--skip-bias\t\t: Don't do bias correction")
  #print("\t--skip-linearity\t: Don't do linearity correction")
  #print("\t--skip-flats\t\t: Don't do flats")
  #print("\t--data-path <opt>\t: NI : Where is the raw data")
  #print("\t--proc-path <opt>\t: NI : Where to write processed data")
  #print("\t--single-chip\t\t: Chip number, default is all.")
  #print("\t--mosaic\t\t: NI : Build Data mosaic")
  #print("\t--clean\t\t\t: NI : Remove files created")
  #print("\n\tNI : Not Implemented (yet)")
  #print("\n\n\tContact: simon@lco.cl\n")
  

if __name__ == '__main__':
  #if "--help" in sys.argv or "-h" in sys.argv:
    #print_help()
    #sys.exit(0)
  if "--data-path" in sys.argv:
    dp_index = sys.argv.index("--data-path")
    new_data_path = sys.argv[dp_index+1]
    if os.path.isdir(new_data_path):
      data_path = new_data_path
    else:
      sys.exit("Invalid Data Path")
  if "--proc-path" in sys.argv:
    pp_index = sys.argv.index("--proc-path")
    new_proc_path = sys.argv[pp_index+1]
    if os.path.isdir(new_proc_path):
      proc_data_path = new_proc_path
    else:
      sys.exit("Invalid Proc Data Path")
  if "--single-chip" in sys.argv:
    index = sys.argv.index("--single-chip")
    if len(sys.argv) > index+1:
      try:
        if int(sys.argv[index+1]) and (0 < int(sys.argv[index+1]) < 5):
          quadrants=[int(sys.argv[index+1])]
          debug("Processing only quadrant %s"%quadrants[0])
          print_spacers("Processing only quadrant %s"%quadrants[0])
        else:
          sys.exit("Argument out of range for chip number \" %s \""%sys.argv[index+1])
      except ValueError:
        sys.exit("Invalid argument for chip number \" %s \""%sys.argv[index+1])
    else:
      sys.exit("Argument missing after --single-chip")
  else:
    print_spacers("Processing all four chips")
  print_spacers("Images")
  image_list	= make_dir_list(data_path)
  debug("%s Images detected in this directory"%len(image_list))
  print_spacers("Classification of files")
  sorted_list	= sort_list(image_list) #bias flat sci
  debug("Processing All, Overscan and Trim.")
  if "--skip-overscan" not in sys.argv:
    print_spacers("Overscan and Trimming")
    overscan_trim(all_list)
  if "--skip-bias" not in sys.argv:
    debug("Producing Master Bias")
    #print bias_list
    print_spacers("Master Bias")
    master_bias	= make_bias(bias_list)
    debug("Applying Bias correction to Flats and Science images")
    print_spacers("Bias Correction")
    bias_correct(flat_sci_list,master_bias)
  if "--shutter" in sys.argv:
    print_spacers("Shutter Correction")
    shutter_corr(flat_sci_list)
  if "--skip-linearity" not in sys.argv:
    print_spacers("Linearity Correction")
    linearity_correction(flat_sci_list)
  if "--skip-flats" not in sys.argv:
    for i in range(len(filter_list)):
      print_spacers("Processing Filter %s"%filter_list[i])
      debug("Processing Filter %s"%filter_list[i])
      debug("Processing Flat")
      debug("Flat List %s %s"%(filter_list[i],len(flat_list[i])))
      make_flats(flat_list[i],filter_list[i])
      #####
      debug("Processing Science")
      debug("Sci List %s %s"%(filter_list[i],len(sci_list[i])))
      flat_correct(sci_list[i],filter_list[i])
  if "--mosaic" in sys.argv:
    for i in range(len(sci_list)):
      if "--single-chip" not in sys.argv:
        mosaic(sci_list[i])
      else:
        debug("Can't do mosaic with single chip")
  if "--clean" in sys.argv:
    to_remove = glob.glob("[t,b,l]_ccd*fits")
    for each in to_remove:
      print(each)
  #print filter_list
  print_spacers("All Done")