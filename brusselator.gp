set terminal pngcairo size 800,600
set output 'brusselator_plot.png'
set title 'Brusselator Model Simulation'
set xlabel 'x'
set ylabel 'y'
set grid
set xrange [0.999949:1.02015]
set yrange [-1:1]
plot 'brusselator_data.txt' using 1:2 every 10 with linespoints title 'Brusselator Data'
