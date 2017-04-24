#set size 1,1
#set terminal postscript eps enhanced 18 size 4.0,3.0
#set lmargin 8.0
#set rmargin 2.5
#set tmargin 3.5
#set bmargin 3.5
set terminal png

datafile = "output.dat"
set xlabel "x"
set ylabel "z"

set title "Density"
set output "density_2D.png"
plot datafile   using 1:2:4 with image 

set title "Pressure"
set output "pressure_2D.png"
plot datafile   using 1:2:5 with image 

set title "Internal Energy"
set output "eint_2D.png"
plot datafile   using 1:2:6 with image 
