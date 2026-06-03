# Build the basal_mixing documentation PDF.
# Usage: `make` (or `make basal_mixing.pdf`); `make clean` to remove aux files.

LATEX     := pdflatex
LATEXOPTS := -interaction=nonstopmode -halt-on-error

.PHONY: all clean watch

all: basal_mixing.pdf

basal_mixing.pdf: basal_mixing.tex
	$(LATEX) $(LATEXOPTS) $<
	$(LATEX) $(LATEXOPTS) $<     # second pass for refs/labels

clean:
	rm -f *.aux *.log *.out *.toc *.synctex.gz

distclean: clean
	rm -f basal_mixing.pdf

# Quick edit loop with latexmk (if available)
watch:
	latexmk -pdf -pvc -interaction=nonstopmode basal_mixing.tex
