EPStoPDF = epstopdf

ALL: UserGuide.ps UserGuide.pdf 

eps_source = $(wildcard *.eps) 

pdf_source = $(eps_source:.eps=.pdf)

section_sources = Introduction/Introduction.tex \
                  GettingStarted/GettingStarted.tex \
	          Inputs/Inputs.tex \
	          UnitsAndConstants/Units.tex \
                  Particles/Particles.tex \

UserGuide.dvi: UserGuide.tex $(eps_source) $(section_sources) $(eps_source)
	latex UserGuide.tex < /dev/null
	latex UserGuide.tex < /dev/null
	latex UserGuide.tex < /dev/null
	latex UserGuide.tex < /dev/null

UserGuide.pdf: UserGuide.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex UserGuide.tex < /dev/null
	pdflatex UserGuide.tex < /dev/null
	pdflatex UserGuide.tex < /dev/null
	pdflatex UserGuide.tex < /dev/null

pdf:	UserGuide.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) UserGuide.dvi
	$(RM) UserGuide.blg UserGuide.log UserGuide.toc
	$(RM) UserGuide.aux UserGuide.ps UserGuide.bbl UserGuide.pdf
	$(RM) *~

.PHONY: clean
