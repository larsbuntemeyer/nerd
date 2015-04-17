set size 1,1
set terminal postscript eps enhanced 18 size 4.0,3.0
set lmargin 8.0
set rmargin 2.5
set tmargin 3.5
set bmargin 3.5

exact_solution_file = "e1rpex.out"
datafile = "output.dat"
unset key
set xlabel "x"
set format y "%2.1f" 

set style line 1 lt 1 lw 1.0 lc 3 
set style line 2 lt 2 lw 1.0

set ylabel "{/Symbol r}"
set yrange [0:1.1]
set title "Density"
set output "sod_shock_density.eps"
plot exact_solution_file using 1:2 w l ls 1, \
     datafile          using 1:2

set ylabel "v"
set yrange [-0.1:1.6]
set title "Velocity"
set output "sod_shock_velocity.eps"
plot exact_solution_file using 1:3 w l ls 1, \
     datafile          using 1:3

set ylabel "p"
set yrange [0:1.1]
set title "Pressure"
set output "sod_shock_pressure.eps"
plot exact_solution_file using 1:4 w l ls 1, \
     datafile          using 1:4

