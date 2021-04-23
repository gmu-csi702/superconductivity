 set title "Band structure of  SbH_3  BCC 6.9      "
  set ylabel "Energy (Ry)" 
 set arrow     1 from   1.000, graph 0 to   1.000   , graph 1 nohead 
 set arrow     2 from   1.707, graph 0 to   1.707   , graph 1 nohead 
 set arrow     3 from   2.414, graph 0 to   2.414   , graph 1 nohead 
 set arrow     4 from   3.280, graph 0 to   3.280   , graph 1 nohead 
 set arrow     5 from   3.780, graph 0 to   3.780   , graph 1 nohead 
 set arrow     6 from   4.280, graph 0 to   4.280   , graph 1 nohead 
 set arrow     7 from   5.146, graph 0 to   5.146   , graph 1 nohead 
set encodi
 set xtics  ("{/Symbol G}" 0.,  "{/Symbol D}"   0.500           ,"H"   1.000           ,"G"   1.354           ,"N"   1.707 ,"{/Symbol S}"   2.061 ,"{/Symbol G}"   2.414 ,"{/Symbol L}"   2.847           ,"P"   3.280           ,"D"   3.530           ,"N"   3.780           ,"D"   3.530           ,"P"   4.280           ,"F"   4.713           ,"H"   5.146  )
 efermi =   1.1035400000000000     
set arrow      8 from 0.000  ,   1.10354 to   5.146  ,   1.10354   nohead
set yrange [  -0.2500000000000000      :   2.4000000000000000      ]
 set xrange [ 0.0    :   5.1462643699419637      ]
 set term postscript eps enhanced solid
 set output "band.eps" 
 plot "banddata" u 1:2 t "" w l lt 1,\
      "banddata" u 1:3 t "" w l lt 1,\
      "banddata" u 1:4 t "" w l lt 1,\
      "banddata" u 1:5 t "" w l lt 1,\
      "banddata" u 1:6 t "" w l lt 1,\
      "banddata" u 1:7 t "" w l lt 1,\
      "banddata" u 1:8 t "" w l lt 1,\
      "banddata" u 1:9 t "" w l lt 1,\
      "banddata" u 1:10 t "" w l lt 1 
 set output 
 
