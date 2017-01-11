#!/bin/bash
#File to build the documenation using sphinx
#Builds if source directory md5sum is different than the saved hash
#Tyler Cope, 2016
#
echo 'Checking all ".py" and ".rst" files for changes...'
echo
#
#Get current and last sum of source directory
FILE_SUM="$(cat lastMd5sum)"
#Get the python and rst hashes
RST_FILE_SUM=`find ./source -name "*.rst*" -exec md5sum {} \; | sort -k 2 | md5sum`
PY_FILE_SUM=`find ../ -name "*.py*" -exec md5sum {} \; | sort -k 2 | md5sum`
#sum the rst source and python hashes
CURRENT_SUM=`echo $RST_FILE_SUM $PY_FILE_SUM | md5sum`
#chop off the end  
CURRENT_SUM=${CURRENT_SUM:0:32}
FILE_SUM=${FILE_SUM:0:32}
#
echo "Current hash = "$CURRENT_SUM
echo "Last hash    = "$FILE_SUM
#
#Write current sum to last sum file
echo -n $CURRENT_SUM > lastMd5sum
#
#If sums dont match build the docs
if [ "$CURRENT_SUM" != "$FILE_SUM" ]
        then
                echo
                echo "File changes found, building docs..."
                echo 
                sphinx-apidoc -f -o source/ ..
                sleep 0.25
                wmctrl -r "Ocelot Doc Builder" -e '0,5000,0,-1,-1'
                make html
                #make latexpdf
                #firefox ./build/html/index.html &
        else
                sleep 0.25
                echo "Source files unchanged, aborting build"
fi
