#!/usr/bin/gnuplot -persist

#set terminal qt size 8192, 8192
set terminal png truecolor size 1024, 1024
#set terminal pdf size 2048, 2048
#unset border
#unset xtics
#unset ytics
#unset ztics
set xlabel "X" font "Linux Libertine O 24"
set ylabel "Y" font "Linux Libertine O 24"
set zlabel "Z" font "Linux Libertine O 24"
set style fill transparent solid 0.20 border
unset key
#do for [AZ=0:360:20] {
    set view 32, 293, 0.9
    set view equal xyz
    outfile = sprintf('landscape.png')
    set output outfile
    splot filename pt 7 ps 0.2 lc rgb 'black'
#}
