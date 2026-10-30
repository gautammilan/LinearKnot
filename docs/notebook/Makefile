# Build the LinearKnot project notebook PDF (pdflatex x2 for refs/TOC).
TEX = linearknot_notebook.tex
PDF = linearknot_notebook.pdf
LATEX = pdflatex
LATEXFLAGS = -interaction=nonstopmode -halt-on-error

.PHONY: all clean

all: $(PDF)

$(PDF): $(TEX)
	$(LATEX) $(LATEXFLAGS) $(TEX)
	$(LATEX) $(LATEXFLAGS) $(TEX)

clean:
	rm -f *.aux *.log *.out *.toc *.pdf
