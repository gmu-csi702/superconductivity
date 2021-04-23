
# Name: dosplot.gnu
# Objective: plot density of states 


EF =system("head -1 ferv_omega | awk '{print $1}'")
EF = EF+0


set title 'Electronic Density of States, a=6.8 a.u.'
set xlabel "Energy(Ry)"
set ylabel "States/Ry/Atom"
set xrange[-1.5:1.5]
set key left top Right noreverse box linetype -2 linewidth 1.000 samplen 4 spacing 1 width 0
set nolabel
set label 1 "{/Symbol e}_F" at 0.0, 10, 0 left norotate
set noarrow
set arrow 1 from 0.0, graph 0, 0 to 0.0, graph 1, 0  nohead linetype 3 linewidth 1.000
#set nolinestyle
unset style line
set nologscale
set offsets 0, 0, 0, 0
set pointsize 1
set encoding default
set label 2 "IH_3" at graph 0.9, 0.9
#Plot total DOS and its angular momentum decompositions
plot "dosdat.bin.plot" u ($2-EF):4 t "Total DOS" w l,     "dosdat.bin.plot" u ($2-EF):5 t "DOS---s" w l,     "dosdat.bin.plot" u ($2-EF):6 t "DOS---p" w l,     "dosdat.bin.plot" u ($2-EF):7 t "DOS---d" w l

set term postscript eps size 10.0,7.0 enhanced color font 'Helvetica,20' linewidth 2
set output "dosdatall.eps"
replot

#Plot total DOS 
plot "dosdat.bin.plot" u ($2-EF):4 t "Total DOS" w l

set label 2 "IH_3" at graph 0.9, 0.9
set arrow 1 from 0.0, graph 0, 0 to 0.0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdattotal.eps"
replot

#Plot total DOS and its angular momentum decompositions
plot "dosdat.bin.plot" u ($2-EF):5 t "DOS---s" w l,     "dosdat.bin.plot" u ($2-EF):6 t "DOS---p" w l,     "dosdat.bin.plot" u ($2-EF):7 t "DOS---eg" w l

set arrow 1 from 0.0, graph 0, 0 to 0.0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdatang.eps"
replot

#Plot S+P angular momentum decompositions
plot "dosdat.bin.plot" u ($2-EF):5 t "DOS---s" w l,     "dosdat.bin.plot" u ($2-EF):6 t "DOS---p" w l

set arrow 1 from 0.0, graph 0, 0 to 0.0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdatsp.eps"
replot

#Plot S angular momentum decompositions
plot "dosdat.bin.plot" u ($2-EF):5 t "DOS---s" w l

set arrow 1 from 0.0, graph 0, 0 to 0.0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdats.eps"
replot

#Plot P angular momentum decompositions
plot "dosdat.bin.plot" u ($2-EF):6 t "DOS---p" w l

set arrow 1 from 0.0, graph 0, 0 to 0.0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdatp.eps"
replot

#Plot D angular momentum decompositions
plot "dosdat.bin.plot" u ($2-EF):7 t "DOS---d" w l

set arrow 1 from 0.0, graph 0, 0 to 0.0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdatd.eps"
replot

set output      


set print "../e_fermi" append
print 6.8, "\t" , EF 
#    EOF
