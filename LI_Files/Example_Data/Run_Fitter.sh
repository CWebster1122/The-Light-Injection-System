# Bash script to run the ROOT fitting macro over the LI data files to
# produce the fit parameters file(s) that will be analyzed by the
# LI_Analysis.C macro. 

temp_file=temp.root
parameter_file_1=Parameters1.dat
parameter_file_2=Parameters2.dat
parameter_file_3=Parameters3.dat

# NOTE: This script as well as the fitter.C macro must be located or
# copied into the directory with the raw data files
filenames=`ls *6` # Can me modified to fit other naming structures

# If the output fit parameters file(s) exist, then remove them to not
# just append to them.
if [ -f $parameter_file_1 ]
then
    rm $parameter_file_1
fi
if [ -f $parameter_file_2 ]
then
    rm $parameter_file_2
fi
if [ -f $parameter_file_3 ]
then
    rm $parameter_file_3
fi

# Loop over all the raw data files and run the fitter on each. In each
# iteration, the current raw data file is copied to temp.root which is
# then fitted. The Fitter.C macro appends the resultant fit parameters
# to a new line in the parameter file(s) and this script then adds the
# time step to the end of each such new line.
date

for name in $filenames; do

    echo " *** Fitting " $name

    # Copy data file to temp file
    if [ -f $temp_file ]
    then
      rm $temp_file
    fi
    cp $name $temp_file

    # Extract the time stamp (in seconds since the start of the year of
    # the first run).
    # time_root=`echo $name | awk -F "data" '{print $2}' | awk -F "_" '{print $2}' | awk -F ".root" '{print $1}'`
    time_root=`echo $name | awk -F ".root" '{print $1}'`

    # Execute the ROOT fitting macro
    root -l -b -q Fitter.C > /dev/null

    # Append the time stamp to the recently added line the parameter file(s)
    echo "$(cat $parameter_file_1) ${time_root}" > $parameter_file_1
    echo "$(cat $parameter_file_2) ${time_root}" > $parameter_file_2
    echo "$(cat $parameter_file_3) ${time_root}" > $parameter_file_3

    # Lastly, the lines corresponding to bad fits are removed
    sed -i '' '/ERROR/d' ./Parameters1.dat
    sed -i '' '/ERROR/d' ./Parameters2.dat
    sed -i '' '/ERROR/d' ./Parameters3.dat
date    
done
