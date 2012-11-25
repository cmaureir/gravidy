TEX = $(wildcard *.tex)
PDF = $(TEX:.tex=.pdf)

all: pdf
pdf:  $(PDF)

%.pdf: %.tex
	@rubber --pdf $<
clean:
	@rubber --clean --pdf $(TEX:.tex=)
x:
	@okular $(PDF) &> /dev/null &
.PHONY: pdf clean

