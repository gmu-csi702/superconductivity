#
# This runs interpot
#

#make sure an element ort.compound was given.
if ($#argv == 0) then
  echo 'you must supply an element ort.compound name'
  exit 1
endif


echo 'removing old links'
rm lapwpot.dat

#make symbolic links to datafiles
echo 'make symbolic links to datafiles'
cp  {$argv[1]}pot.lapw lapwpot.dat
if ($status != 0) then
  echo 'error exit in interpot.com - linking in file'
  exit 1
endif
#run interpot
echo 'run interpot'
./interpot
if ($status != 0) then
  echo 'error exit in interpot.com - interpot program failed'
  exit 1
endif
#rename output files
echo 'rename output files'
mv -f apwpot.dat {$argv[1]}.pot

#remove symbolic links
echo 'remove symbolic links'
rm lapwpot.dat
