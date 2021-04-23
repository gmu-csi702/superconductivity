#! /bin/sh 
rm BRH* 
 ~/lapw/lapw4at  > /dev/null
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



