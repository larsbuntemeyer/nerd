#set size 1,1
#set terminal postscript eps enhanced 18 size 4.0,3.0
#set lmargin 8.0
#set rmargin 2.5
#set tmargin 3.5
#set bmargin 3.5
set terminal png

datafile = "output.dat"
#unset key
set xlabel "x"
set ylabel "z"

#set style line 1 lt 1 lw 1.0 lc 3 
#set style line 2 lt 2 lw 1.0

#set pm3d map

#set view map

#set ylabel "{/Symbol r}"
set title "Density"
set output "sod_shock_density_2D.png"
splot datafile   using 1:2:4 

