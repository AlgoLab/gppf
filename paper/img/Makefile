sources=$(wildcard *.tex)
pdfs=$(sources:.tex=.pdf)
pdf: $(pdfs)
	

%.pdf : %.tex 
	pdflatex -shell-escape $<


clean :
	latexmk -c && rm -f $(pdfs)

