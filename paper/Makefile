pdfs= single-cell-ilp.pdf
bibs=$(wildcard *.bib)

all : $(pdfs)

%.pdf : %.tex $(bibs) figs
	latexmk -pdf $<

figs :
	cd img && make

clean :
	latexmk -c && rm -f $(pdfs) && cd img && make clean
