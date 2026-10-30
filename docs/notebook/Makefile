# Host texlive-cm TFMs may be missing; set TEXMFHOME to a tree with
# fonts/tfm/public/cm/cmr10.tfm (e.g. $(HOME)/texmf after seeding CM TFMs).
TEXMFHOME ?= $(HOME)/texmf
export TEXMFHOME

.PHONY: all clean
all: linearknot_notebook.pdf

linearknot_notebook.pdf: linearknot_notebook.tex
	pdflatex -interaction=nonstopmode $<
	pdflatex -interaction=nonstopmode $<

clean:
	rm -f *.aux *.log *.out *.toc missfont.log
