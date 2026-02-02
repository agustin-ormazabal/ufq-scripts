reset

set terminal pngcairo size 1800,700
set output "perfiles_wt.png"

set lmargin 11 
#set rmargin 8
set multiplot layout 1,3

set size square
set xrange[12:31.5]
set yrange[0:80]

set grid ytics lw 1

set key at 30.5,12 font "Verdana,15" tc rgb "#555555"
set ytics font "verdana,20" tc rgb "#555555"
set xtics font "verdana,15" tc rgb "#555555"

#############
# WILD TYPE #
#############

unset title
unset ylabel
set origin 0.1,0.047
set style fill transparent solid 0.3 noborder

p "hcnA_2.out" u ($1-1.2):($2+0.8) smooth csplines lw 5.5 lc rgb "#332d63" t "RsmE-(SL0)_{2}", "hcnA_2_std.out" u ($1-1.2):($2+0.8+$3):($2+0.8-$3) with filledcurves lc rgb "#332d63" notitle, "hcnA_1_std.out" u 1:($2+0.5) smooth csplines lw 6 lc rgb "#c87d7d" t "RsmE-(SL0)_{1}", "hcnA_1_std.out" u 1:($2+0.5+$3):($2+0.5-$3) with filledcurves lc rgb "#c87d7d" notitle


set xrange[12.28:30.8]

set origin 0.38333,0.047
set key at 30.5,12 font "Verdana,15" tc rgb "#555555"
set ytics font "verdana,0" tc rgb "#FFFFFF"
set grid ytics lw 1

p "SL2_2_std.out" u ($1+0.5):($2+0.5) smooth csplines lw 6 lc rgb "#332d63" t "RsmE-(SL2)_{2}", "SL2_2_std.out" u ($1+0.5):($2+0.5-$3):($2+0.5+$3) with filledcurves lc rgb "#332d63" notitle, "SL2_1.out" u ($1+0.5):($2+0.1) smooth csplines lw 6 lc rgb "#c87d7d" t "RsmE-(SL2)_{1}", "SL2_1_std.out" u ($1+0.5):($2+0.1+$3):($2+0.1-$3) with filledcurves lc rgb "#c87d7d" notitle

set xrange[12.28:25.5]
set origin 0.6666,0.047
set key at 24.5,69.25 font "Verdana,15" tc rgb "#555555"
set ylabel "PMF (kcal/mol)" font "Verdana,30" tc rgb "#555555" offset -110,0,0
set label "Reaction coordinate (Å)" font "Verdana,30" tc rgb "#555555" at -1.5,-14
set grid ytics lw 1 

p "SL3_2_std.out" u ($1-0.5):2 smooth csplines lw 6 lc rgb "#332d63" t "RsmE-(SL3)_{2}", "SL3_2_std.out" u ($1-0.5):($2-$3):($2+$3) with filledcurves lc rgb "#332d63" notitle, "SL3_1.out" u 1:2 smooth csplines lw 6 lc rgb "#c87d7d" t "RsmE-(SL3)_{1}", "SL3_1_std.out" u 1:($2+$3):($2-$3) with filledcurves lc rgb "#c87d7d" notitle

unset multiplot
