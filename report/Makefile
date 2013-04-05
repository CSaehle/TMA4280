LATEX = pdflatex
BIBTEX = bibtex
LATEXFLAGS = -shell-escape
TEMPS = toc bbl blg log out lof lot lol aux pyg
REPORT = tma4280_report

.PHONY: clean purge remake
.SUFFIXES: .tex .pdf

.tex.pdf:
	$(LATEX) $(LATEXFLAGS) $<
#	$(BIBTEX) $*.aux
#	$(LATEX) $(LATEXFLAGS) $<
	$(LATEX) $(LATEXFLAGS) $<

all: | purge $(REPORT).pdf

clean:
	for s in $(TEMPS); do rm -f $(REPORT).$$s; done

purge: clean
	-rm -f $(REPORT).pdf
