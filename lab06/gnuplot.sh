set term post colour solid enhanced  # format rysunku
set term png size 1024,768 font "Times Roman,19"

set output "funkcje.png"
#set title "Wychylenie x(t)"
set xlabel "x"
set ylabel "f(x)"
#set key center center
set yrange [-1 : 2.5]
set xrange [0.9 : 3.7]

set grid x y				 # siatka pomocnicza
set style line 1 lt 1 lw 1 linecolor rgb '#4169e1' pt 7 ps 1.5
set style line 2 lt 1 lw 1 linecolor rgb '#8b0000' pt 7 ps 1.5

plot 'funkcje.txt' using 1:2 ls 1 w l title 'f(x)', 'funkcje.txt' using 1:3 ls 2 w l title '{f}^{\''}(x)'

# w p   == with points
# w l   == with line
# t "dt = 0.1" == title "dt = 0.1"
