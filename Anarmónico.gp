clear
set xl "x"
set yl "{/Symbol y}"
set grid
p "anarmonico.dat" u 1:2 w l ls 6 t "{/Symbol y}_L", "anarmonico.dat" u 3:4 w l ls 7 t "{/Symbol y}_R"
