#!/bin/bash

# perl ./parselog.pl result.build.txt |tee result.dat
perl ./parselog.pl ${1:-result.build_no_elimination.txt} |tee result.dat

#gnuplot <<EOF
gnuplot -persist <<EOF
set title "Simulated Watershed Run"
set output "result.eps"
 plot "result.dat" index 0 using 2:3 w linespoints title "PQ NOP", \
 "result.dat" index 1 using 2:3 w linespoints title "PQ ConversionScanOrderIndexToCoordinate", \
  "result.dat" index 2 using 2:3 w linespoints title "PQ SOI on queue", \
 "result.dat" index 3 using 2:3 w linespoints title "Turbo NOP", \
 "result.dat" index 4 using 2:3 w linespoints title "Turbo ConversionScanOrderIndexToCoordinate", \
  "result.dat" index 5 using 2:3 w linespoints title "Turbo SOI on queue"
set ylabel "time [ms]"
set xlabel "payload [bytes]"
 pause -1
# replot
EOF
# gzip result.eps
#  plot "result.dat" index 0 using 3:xticlabel(2) w linespoints title "NOP", \
# "result.dat" index 1 using 3:xticlabel(2) w linespoints title "ConversionScanOrderIndexToCoordinate", \
#  "result.dat" index 2 using 3:xticlabel(2) w linespoints title "ConversionScanOrderIndexToCoordinateAndBack"

