#! /bin/sh 
rm BRH* 
 ~/lapw/ih3/lat8.6/lapw4at-8.6  > /dev/null
if [ -f "ATSPDT" ]
then
	rm ATSPDT
fi
if [ -f "fort.16" ]
then
	rm fort.16
fi
if [ -f "AFOP" ]
then
	rm AFOP
fi
if [ -f "CDN0" ]
then
	rm CDN0
fi
rm RECIP
rm CORECG
rm EIG*
rm IO*
rm PLFILE
rm W*
rm fort*



