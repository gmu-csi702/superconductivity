#! /bin/sh 
rm BRH* 
 ~/lapw/ih3/lat6.1/lapw4at-6.1  > /dev/null
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



