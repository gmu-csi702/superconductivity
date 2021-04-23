#!/bin/bash 
el="SbH_3"
pt=18
rm e_fermi

for i in $(seq 7.0 0.2 8.4)
do
mkdir -p "lat"$i
cd lat$i

cp ~/lapw/sbh3/bcc/lat_src/* . 
cp INFILE.iter INFILE
sed -i "5 c\ $i " INFILE

rm BRH* 
~/lapw/lapw4at > /dev/null
 
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



cdn1=$(ls -t CDN* | head -1)

mv CDN1 CDN1.iter
mv $cdn1 CDN1

mv INFILE INFILE.iter
mv INFILE.run INFILE
sed -i "5 c\ $i " INFILE
sed -i "39 c\ -1.500     2.000" INFILE

rm BRH* 
~/lapw/lapw4at > /dev/null 
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



cp ../dos_src/{qlmtconvert,dosapwn} .

./qlmtconvert

cat > dosdat.out << EOF
   11    4  101
Sb  H3                                      2
  8.00 54.00
   1   1
   $i         1    2
    8    9   11    1    0
  -1.50000   0.00200   2.00000   0.00000   0.00000   0.00000   0.00000

EOF

./dosapwn


cat > dosplot.gnu << EOF

# Name: dosplot.gnu
# Objective: plot density of states 


EF =system("head -1 ferv_omega | awk '{print \$1}'")
EF = EF+0


set title "Electronic Density of States"
set xlabel "Energy(Ry)"
set ylabel "States/Ry/Atom"

set key left top Right noreverse box linetype -2 linewidth 1.000 samplen 4 spacing 1 width 0
set nolabel
set label 1 "{/Symbol e}_F" at EF, $pt, 0 left norotate
set noarrow
set arrow 1 from EF, graph 0, 0 to EF, graph 1, 0  nohead linetype 3 linewidth 1.000
#set nolinestyle
unset style line
set nologscale
set offsets 0, 0, 0, 0
set pointsize 1
set encoding default
set label 2 "$el" at graph 0.9, 0.9
#Plot total DOS and its angular momentum decompositions
plot "dosdat.bin.plot" u 2:4 t "Total DOS" w l,\
     "dosdat.bin.plot" u 2:5 t "DOS---s" w l,\
     "dosdat.bin.plot" u 2:6 t "DOS---p" w l,\
     "dosdat.bin.plot" u 2:7 t "DOS---d" w l

set term postscript eps size 10.0,7.0 enhanced color font 'Helvetica,20' linewidth 2
set output "dosdatall.eps"
replot

#Plot total DOS 
plot "dosdat.bin.plot" u 2:4 t "Total DOS" w l

set label 2 "$el" at graph 0.9, 0.9
set arrow 1 from EF, graph 0, 0 to EF, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdattotal.eps"
replot

#Plot total DOS and its angular momentum decompositions
plot "dosdat.bin.plot" u 2:5 t "DOS---s" w l,\
     "dosdat.bin.plot" u 2:6 t "DOS---p" w l,\
     "dosdat.bin.plot" u 2:7 t "DOS---eg" w l

set arrow 1 from EF, graph 0, 0 to EF, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdatang.eps"
replot

#Plot S+P angular momentum decompositions
plot "dosdat.bin.plot" u 2:5 t "DOS---s" w l,\
     "dosdat.bin.plot" u 2:6 t "DOS---p" w l

set arrow 1 from EF, graph 0, 0 to EF, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdatsp.eps"
replot

#Plot S angular momentum decompositions
plot "dosdat.bin.plot" u 2:5 t "DOS---s" w l

set arrow 1 from EF, graph 0, 0 to EF, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdats.eps"
replot

#Plot P angular momentum decompositions
plot "dosdat.bin.plot" u 2:6 t "DOS---p" w l

set arrow 1 from EF, graph 0, 0 to EF, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdatp.eps"
replot

#Plot D angular momentum decompositions
plot "dosdat.bin.plot" u 2:7 t "DOS---d" w l

set arrow 1 from EF, graph 0, 0 to EF, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdatd.eps"
replot

set output      


set print "../e_fermi" append
print $i, "\t" , EF 
#    EOF
EOF

gnuplot dosplot.gnu



cp -r ~/lapw/sbh3/bcc/gyo_src/ gyo
cd gyo/

cat > input.dat << EOF
SbH3                 LAPW                   2     
Sb
  435   0.0
  480   10   32  2.2000     10000
H
  435   0.0
  400   10   32  1.2000     10000

EOF

cp ../POT lapwpot.dat

./interpot

cp apwpot.dat pot.in

head -n2 ../dosapw.itp > dos.in
sed -i "3 c\   121.75    500.     0.10    1334.6" dos.in 
tail -n1 ../dosapw.itp >>  dos.in
sed -i "5 c\   1.0079   1750.     0.10    1334.6" dos.in

./gyond 


rm  fort.*

cd ../

cd ../

pt=`expr $pt + 2`
 
done



