set termoption enhanced
set termoption font 'Helvetica,16'
set xlabel '{/Symbol b}'
set ylabel 'E_{P}'
plot 'plaquette.txt' with errorbars title '16^{4} lattice' pointtype 7
