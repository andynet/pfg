.ONESHELL:

.PHONY: build
build:
	latexmk -pdf -shell-escape main.tex

.PHONY: clean
clean:
	rm -r *abs *aux *bbl *blg *fdb_latexmk *fls *log *out *xmpdata _minted*

.PHONY: purge
purge:
	rm main.pdf

