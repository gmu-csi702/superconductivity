set term pdf
set output "sbh3_eta.pdf"

set title "SbH_3 Hopfield Parameter, {/Symbol h}"
set ylabel "{/Symbol h}"
set xlabel "a (A)"
set yrange[0.5:2.5]
plot "hopfield_par.corr" u 1:2 w p pt 5 ps 0.65 t "{/Symbol h}_{Sb}",     "hopfield_par.corr" u 1:($3*3) w p pt 3 ps 0.65 t "{/Symbol h}_H"


