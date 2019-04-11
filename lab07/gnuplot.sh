set term post colour solid enhanced  # format rysunku
set term png size 1024,768 font "Times Roman,19"

set output "n20.png"
#set title "Wychylenie x(t)"
set xlabel "x"
set ylabel "f(x)"
#set key center center
set xrange [-5 : 5]

set grid x y				 # siatka pomocnicza
set style line 1 lt 1 lw 1 linecolor rgb '#4169e1' pt 7 ps 1.5
set style line 2 lt 1 lw 1 linecolor rgb '#8b0000' pt 7 ps 1.5

plot 'n20.dat' using 1:2 ls 1 w l title 'W_{20}(x)', 'exp.dat' using 1:2 ls 2 w l title 'f(x) = exp(-x^2)'

# w p   == with points
# w l   == with line
# t "dt = 0.1" == title "dt = 0.1"
