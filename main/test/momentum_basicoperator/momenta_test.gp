set title "A80 p = 1 test (2 random vectors)"
set term postscript
set out "momenta_test.ps"
set xrange [0:70]
set xlabel "t/a"
set ylabel "Re(C2(p))
plot "./plot.txt" u 1:2 t "(-1,  0,  0)", "./plot.txt" u 1:3 t "( 0, -1,  0)", "./plot.txt" u 1:4 t "( 0,  0, -1)", "./plot.txt" u 1:5 t "( 0,  0,  0)", "./plot.txt" u 1:6 t "( 0,  0,  1)", "./plot.txt" u 1:7 t "( 0,  1,  0)", "./plot.txt" u 1:8 t "( 1,  0,  0)"
