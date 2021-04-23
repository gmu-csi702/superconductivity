#!/bin/bash


#loop over the lattice constants of interest
#for i in $(seq 6.0 0.2 6.8)
for i in  6.6 
do
    #cp -r a_ini  lat$i
    cd lat$i/c11-c12/
    #cd lat$i/c44/

    #for j in e0000  e0005  e001  e002  e003  e004  e005  e006  #e007  e008
    for j in e0000  e0005  e001  e002  e004  e006  
    do
        cd $j
	#pwd
        grep "TOTAL ENERGY" INFO | tail -1 
	cd ../
    done
    cd ../../
done
