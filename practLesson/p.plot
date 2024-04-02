set title "Merry-go-round"
set xlabel "t"
set ylabel "x"
plot "< ./prog2 1 100 0.01" using 1:2 w lines t ""
plot "< ./prog2 1 100 0.5" using 2:3 w lines t ""
pause -1