set term eps
set output "sbh3_pressure.eps"
set title "SbH_3 Pressure"
set xlabel "a (A)"
set ylabel "pressure (Mbar)"
set grid

plot 'pressfit.out' u 1:3 w l  notitle

