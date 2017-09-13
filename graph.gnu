set autoscale
unset log
unset label
set xtic auto
#set xtic 3
#set ytic 1
#set y2tic auto

#set logscale y

#set title "Title"
#set xlabel "Simulation time (min)"
set ylabel "Correlation"

set xlabel font ",22"
set ylabel font ",22"


#set xrange [0 : 15]
#set yrange [0 : 100]


#Legend
set key bottom right
#set key top center
set key font ",22"
set key spacing 1.75
set key width 12
#set key outside

set terminal postscript eps enhanced color font "Times"


set output "correlation-lambda3.eps"
set xlabel "Sensors number"
plot \
"stats/algo1corr_OK_l3.data" using 1:2 with linespoints t "Algo1", \
"stats/optcorr_OK_l3.data" using 1:2 with lines t "Opt", \
"stats/randcorr_OK_l3.data" using 1:2 with lines t "Random"


set output "correlation-lambda4.eps"
set xlabel "Sensors number"
plot \
"stats/algo1corr_OK_l4.data" using 1:2 with linespoints t "Algo1", \
"stats/optcorr_OK_l4.data" using 1:2 with lines t "Opt", \
"stats/randcorr_OK_l4.data" using 1:2 with lines t "Random"


set output "correlation-lambda5.eps"
set xlabel "Sensors number"
plot \
"stats/algo1corr_OK_l4.data" using 1:2 with linespoints t "Algo1", \
"stats/optcorr_OK_l4.data" using 1:2 with lines t "Opt", \
"stats/randcorr_OK_l4.data" using 1:2 with lines t "Random"



























set xlabel "Lambda"
set ylabel "Sensors"
set zlabel "Corr"
#set ticslevel 0.0
#set cntrparam levels auto 1
#set cntrparam levels disc 40,50,60,70,80,90
#set contour case
#set hidden3d back offset 1 trianglepattern 3 undefined 1 altdiagonal bentover
set border 31+32+64+256+512 lw .5
set pm3d
set hidden3d
#set dgrid3d 40,40 qnorm 2
#set dgrid3d 80,80 splines
#set dgrid3d 100,100 hann
#gauss | cauchy | exp | box | hann
#set dgrid3d 140,140 gauss
#set cntrparam bspline
#set pm3d implicit at bstbst
#set pm3d implicit at b

set output "correlation-Algo1-3d.eps"
splot "stats/algo1corr_OK_3D.data" using 1:2:3 title '' with lines

set output "correlation-Opt-3d.eps"
splot "stats/optcorr_OK_3D.data" using 1:2:3 title '' with lines

set output "correlation-Rand-3d.eps"
splot "stats/randcorr_OK_3D.data" using 1:2:3 title '' with lines





















