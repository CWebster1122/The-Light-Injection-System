#!/bin/bash
 n_of_runs=$1
 output_name=acquisition.root
 new_name_prefix=/home/jenny/CAMAC_DAQ/camac_benton/LI_Data/Chris_Data/new_reptest_Bi/

 
 echo " I will run over " $n_of_runs " runs"

 I=$n_of_runs
 #for(( i=0, rep=230; i<I, rep<1100; i++, rep+=200)); do
 for(( i=0; i<I; i++)); do
  echo "doing run " $i
  month=`date | awk '{print $2}'`
  day=`date | awk '{print $3}'`
  hour=`date | awk '{print $4}'`

  echo month $month day $day hour $hour 
  initial_time=`date +%s`
  
  #cd /home/jenny/pulser/pulse/pulse 
  #./pulse -s
  #./pulse -c -h 168 -r $rep -w 1                                      
  #./pulse -c -h $ph
  #cd /home/jenny/CAMAC_DAQ/camac_benton 
  ./camacDaq  
  new_full_name=$new_name_prefix$initial_time.root

  echo file name $new_full_name
  mv $output_name $new_full_name
 done

 echo all done
