#!/bin/sh

rm *.aux *.bbl *.blg *.log *.out *.toc
pdflatex cubfits-guide.Rnw
bibtex cubfits-guide
pdflatex cubfits-guide.Rnw
pdflatex cubfits-guide.Rnw
pdflatex cubfits-guide.Rnw
rm *.aux *.bbl *.blg *.log *.out *.toc

mv -f *.pdf ../inst/doc/
cp -f *.Rnw ../inst/doc/
