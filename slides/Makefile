slides.pdf: slides.tex $(wildcard graphics/*.pdf) references.bib samuel.tex niklas.tex
	latexmk -pdf slides.tex

.PHONY: clean
clean:
	rm -f *.aux *.log *.bbl slides.pdf *.nav *.out *.snm *.toc *.fdb_* *.fls *.blg
