set term eps
set output "sbh3_etot.eps"
set title "SbH_3 Total Energy"
set xlabel "a (A)"
set ylabel "energy (Ry)"

plot 'energyfit.out' u 1:3 w l  notitle

