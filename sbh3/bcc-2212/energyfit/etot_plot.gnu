set term eps
set output "sbh3_etot_lattice.eps"
set title "SbH_3 Total Energy"
set xlabel "a (A)"
set ylabel "Energy/Atom (Ry)"
set xrange[6.5:8.5]
set xtics 6.4,0.2,8.4
plot 'etot.sbh3' u 1:2 w p pt 4 t "Calculated Energies", "energyfit.out" u 1:3 w l t "Fit Energies"

