set term post colour solid enhanced # format rysunku
set term png size 1024,768 font "Times Roman,19"

set output "k8.png"
set title "k = 8, N = 2^k próbek wejściowych"
set xlabel "t"
set ylabel "f(t)"
#set key center center
#set xrange [-5 : 5]

set grid x y				 # siatka pomocnicza
set style line 1 lt 1 lw 3 linecolor rgb '#4169e1' pt 7 ps 1.5
set style line 2 lt 1 lw 1 linecolor rgb '#8b0000' pt 7 ps 1.5
set style line 3 lt 1 lw 3 linecolor rgb '#555555' pt 7 ps 1.5
set style line 4 lt 1 lw 1 linecolor rgb '#ae2100' pt 7 ps 1.5

plot 'fdelta.dat' using 1:2 ls 4 w l title 'f(t_i) = f_0(t_i) + \Delta', 'output.dat' using 1:2 ls 1 w l title 'f(t_i)*g(t_i)', 'f.dat' using 1:2 ls 3 w l title 'f_0(t)'

# w p   == with points
# w l   == with line
# t "dt = 0.1" == title "dt = 0.1"
