gnuplot plot.gpl
epstopdf plot-inc.eps
pdflatex plot.tex
open plot.pdf
pdftocairo -png plot.pdf
