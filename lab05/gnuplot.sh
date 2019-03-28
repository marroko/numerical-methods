set term post colour solid enhanced  # format rysunku
set term png size 1024,768 font "Times Roman,19"

set output "lambdas.png"
#set title "Wychylenie x(t)"
set xlabel "Numer iteracji"
set ylabel "{/Symbol l}_{k}"
set yrange [0 : 0.6]
set encoding utf8

set grid x y				 # siatka pomocnicza
set style line 1 lt 1 lw 1 linecolor rgb '#4169e1' pt 7 ps 1.5
set style line 2 lt 1 lw 1 linecolor rgb '#8b0000' pt 7 ps 1.5
set style line 3 lt 1 lw 1 linecolor rgb '#8bffff' pt 7 ps 1.5
set style line 4 lt 1 lw 1 linecolor rgb '#abcdef' pt 7 ps 1.5
set style line 5 lt 1 lw 1 linecolor rgb '#ffbaf5' pt 7 ps 1.5
set style line 6 lt 1 lw 1 linecolor rgb '#459ffa' pt 7 ps 1.5
set style line 7 lt 1 lw 1 linecolor rgb '#69cc21' pt 7 ps 1.5

plot 'dane.dat' using 1:3 ls 2 w lp title '{/Symbol l}_{2}','dane.dat' using 1:4 ls 3 w lp title '{/Symbol l}_{3}','dane.dat' using 1:5 ls 4 w lp title '{/Symbol l}_{4}','dane.dat' using 1:6 ls 5 w lp title '{/Symbol l}_{5}','dane.dat' using 1:7 ls 6 w lp title '{/Symbol l}_{6}','dane.dat' using 1:8 ls 7 w lp title '{/Symbol l}_{7}'

# w p   == with points
# w l   == with line
# t "dt = 0.1" == title "dt = 0.1"
