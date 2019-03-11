set term post colour solid  # format rysunku
set term png size 1024,768 font "Times Roman,19"

set output "wychylenie.png"
#set title "Wychylenie x(t)"
set xlabel "t"
set ylabel "x(t)"
set grid 					 # siatka pomocnicza
set style line 1 lt 1 lw 1 linecolor rgb '#4169e1' pt 7 ps 0.5
set style line 2 lt 1 lw 1 linecolor rgb '#8b0000'

plot "dane.dat" ls 1 w p title "x(t), dt = 0.1", cos(x) ls 2 w l title "cos(t)"

# w p   == with points
# w l   == with line
# t "dt = 0.1" == title "dt = 0.1"
