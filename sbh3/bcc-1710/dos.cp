#!/bin/bash
mkdir -p dos_files
for i in $(seq 5.6 0.1 8.4)
do
  if [ -e   "lat"$i"/dos/dosdatall.eps" ]; then
    cp  "lat"$i"/dos/dosdatall.eps" dos_files/dos.$i.eps
  fi
done  

