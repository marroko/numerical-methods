set term post colour solid  # format rysunku
set term png size 1024,768 font "Times Roman,19"

set output "double.png"
#set title "Wychylenie x(t)"
set xlabel "k (numer iteracji)"
set ylabel "||r_{k}||_{2}"
set y2label "||x_{k}||_{2}"
unset log y2
set logscale y
set y2range [4695 : 4725]
set ytics nomirror
set y2tics
set tics out
set grid x y				 # siatka pomocnicza
set key right center
set style line 1 lt 1 lw 1 linecolor rgb '#4169e1' pt 7 ps 0.5
set style line 2 lt 1 lw 1 linecolor rgb '#8b0000' pt 7 ps 0.5
set style line 3 lt 1 lw 2 linecolor rgb '#abcdef' pt 7 ps 0.5
set style line 4 lt 1 lw 2 linecolor rgb '#fedcba' pt 7 ps 0.5

plot 'double.dat' using 1:2 ls 3 w l title '||r_{k}||_{2}, x_{0} = 1' axes x1y1, 'double.dat' using 1:4 ls 4 w l title '||x_{k}||_{2}, x_{0} = 1' axes x1y2, 'double.dat' using 1:2 ls 1 w p title '||r_{k}||_{2}, x_{0} = 0' axes x1y1, 'double.dat' using 1:4 ls 2 w p title '||x_{k}||_{2}, x_{0} = 0' axes x1y2,

# w p   == with points
# w l   == with line
# t "dt = 0.1" == title "dt = 0.1"
