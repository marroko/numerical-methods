set term post colour solid enhanced # format rysunku
set term png size 1024,768 font "Times Roman,19"

set output "n101.png"
set title "N = 101"
set xlabel "x"
set ylabel "g(x)"
#set key center center
#set xrange [-5 : 5]

set grid x y				 # siatka pomocnicza
set style line 1 lt 1 lw 1 linecolor rgb '#4169e1' pt 7 ps 1.5
set style line 2 lt 1 lw 1 linecolor rgb '#8b0000' pt 7 ps 1.5

plot 'n101.dat' using 1:2 ls 2 w l title 'G(x); N = 101; alfa = 0.5', 'wn101.dat' using 1:2 ls 2 w p title 'g(x_j)'

# w p   == with points
# w l   == with line
# t "dt = 0.1" == title "dt = 0.1"
