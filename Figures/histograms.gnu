reset

set terminal pngcairo size 1800,700
set output "histo_wt.png"

set lmargin 8 
#set rmargin 8
set multiplot layout 1,3

set size square

set grid xtics lw 1

set yrange[0:0.032]
set key at 425,0.036 font "Verdana,13" tc rgb "#555555"
set ytics font "verdana,16" tc rgb "#555555"
set xtics 20 font "verdana,14" tc rgb "#555555"   

#############
# WILD TYPE #
#############

unset title
unset ylabel
set origin 0.1,0.0
set style fill transparent solid 0.3 noborder
set xrange[221:419]

p "histo_hcna_2.dat" u 1:3 smooth csplines  lw 5.5 lc rgb "#332d63" t "RsmE-(SL0)_{2}", "histo_hcna_1.dat" u 1:3 smooth csplines lw 6 lc rgb "#c87d7d" t "RsmE-(SL0)_{1}"

set origin 0.38333,0.0
set key at 405,0.036 font "Verdana,13" tc rgb "#555555"
set ytics font "verdana,0" tc rgb "#FFFFFF"
set label "Number of contacts" font "Verdana,30" tc rgb "#555555" at 220,-0.005
set xrange[201:399]

p "histo_sl2_2.dat" u 1:3 smooth csplines lw 6 lc rgb "#332d63" t "RsmE-(SL2)_{2}", "histo_sl2_1.dat" u 1:3 smooth csplines lw 6 lc rgb "#c87d7d" t "RsmE-(SL2)_{1}"

unset label
set origin 0.6666,0.0
set key at 325,0.036 font "Verdana,15" tc rgb "#555555"
set ylabel "Probability" font "Verdana,30" tc rgb "#555555" offset -105,0,0

set xrange[121:319]
p "histo_sl3_2.dat" u 1:3 smooth csplines lw 6 lc rgb "#332d63" t "RsmE-(SL3)_{2}", "histo_sl3_1.dat"  u 1:3 smooth csplines lw 6 lc rgb "#c87d7d" t "RsmE-(SL3)_{1}"


unset multiplot
