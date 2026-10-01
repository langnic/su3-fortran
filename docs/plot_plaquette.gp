set termoption enhanced
set termoption font 'Helvetica,16'
set xlabel '{/Symbol b}'
set ylabel 'E_{P}'
plot '../data/plaquette_beta_scan_16x16.txt' with errorbars title '16^{4} lattice' pointtype 7
