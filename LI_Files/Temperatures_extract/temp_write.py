#****************************************************************************************************************************************************
# This code is meant to process the temperature data for the LI system, and output the corresponding data to the file "Temp_Output". 
# "Temp_Output" is then read-in to the LI_Analysis.C macro.
# The logic of this is heavily reliant on the excel structure in which the data was recorded. 
# So if the structure changes, this code will break pretty easily. 
# However, it is generaly an easy/quick fix if you know what the format of the excel file changed to.
# P.S. This file was originally formatted in such a way as to match the incorrect time stamps issued by the original doit_w_dates.sh. 
# It has since been changed to write the correct time stamps to "Temp_Output"
#
# * This specific file is for use with one temperature data file ("tempya.txt" in this case) *
#*****************************************************************************************************************************************************

#!/usr/bin/env python
import sys
import os
import time as t
#start_time = t.mktime((2016,1,1,0,0,0,0,1,-1)) #time stamp of January 1, 2016    # This is for forcing the temp timestamp to match up with incorrect data timestamp
#print("starting time: " + str(t.localtime(start_time)))    # for temp timestamp to match incorrect data timestamps
if os.path.exists("Temp_Output"):
  os.remove("Temp_Output")
  output_file = open("Temp_Output", "a")
else:
  print("The file does not exist")
  output_file = open("Temp_Output", "a")

time_array = []
original_temp_data = open("tempya.txt", 'rU')
line_number = 1
for line in original_temp_data:
  if(line_number == 1):
    line_number = line_number + 1 #this skips the first row which are just column labels
    continue
  line_split = list(map(str.strip, line.split("\t")))
  #print line_split
  print((line_split[1])+(line_split[2]))
  file_time = t.strptime(line_split[1]+line_split[2],'%m/%d/%y%H:%M:%S')
  #print(file_time)
  the_time = t.mktime(file_time)
  #time_elapsed = t.mktime(file_time) - start_time    # for temp timestamp to match incorrect data timestamps
  #print(str(t.mktime(time_elapsed)) + " " + str((t.mktime(file_time) - start_time)))    # for temp timestamp to match incorrect data timestamps
  #print(str(t.mktime(file_time)) + " " + str((t.mktime(file_time) - start_time)))    # for temp timestamp to match incorrect data timestamps
  time_array.append(the_time)

  #print("Readout: " + str(t.localtime(the_time)))    # for temp timestamp to match incorrect data timestamps
  #print("Readout: " + str(t.localtime(time_elapsed)))    # for temp timestamp to match incorrect data timestamps
  #print("Readout without subtraction: " + str(t.localtime(t.mktime(file_time))))    # for temp timestamp to match incorrect data timestamps
  line_number = line_number + 1
  #print("\n \n")
  output_file.write(str(the_time) + " " + str(line_split[3]) + " " + str(line_split[5]) + " " + str(line_split[7]) + " " + str(line_split[9]) + "\n")
  #output_file.write(str(time_elapsed) + " " + str(line_split[3]) + " " + str(line_split[5]) + " " + str(line_split[7]) + " " + str(line_split[9]) + "\n"    # for temp timestamp to match incorrect data timestamps)
