PROC=pdflatex
OUTDIR=output
ARGS=--output-directory $(OUTDIR)
FILENAME=presentation

makepdf:
	$(PROC) $(ARGS) $(FILENAME).tex
	cp $(OUTDIR)/$(FILENAME).pdf .
