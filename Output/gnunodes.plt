#!/usr/bin/gnuplot -persist

#set terminal qt size 8192, 8192
set terminal png size 1024, 1024
unset border
unset xtics
unset ytics
unset ztics
unset key
do for [AZ=0:360:20] {
    set view 75, AZ, 1.2
    set view equal xyz
    outfile = sprintf('landscape%d.png',AZ)
    set output outfile
    splot filename pt 7 ps 0.5 lc rgb 'black'
}
