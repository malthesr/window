OUTDIR = out
TMPDIR = build
SECDIR = sections
SUPDIR = supplementary
RDIR = scripts
TIKZDIR = tikz
FIGDIR = figures
TABDIR = tables

MAIN := $(OUTDIR)/main.pdf
SUPPLEMENTARY := $(OUTDIR)/supplementary.pdf
BOTH := $(OUTDIR)/main_supplementary.pdf
DEPS = imports.sty macros.tex

SECNAMES = abstract introduction methods results discussion
SECS := $(addprefix $(SECDIR)/, $(addsuffix .tex, $(SECNAMES)))

FIGNAMES = window.pdf 1d.png 2d.png dadi.png bench_impala.png sim_loglik.png
FIGS := $(addprefix $(FIGDIR)/, $(FIGNAMES))
FIGDEPS := $(wildcard $(RDIR)/plotting/*.R)

TABS := $(TABDIR)/data_train.tex

SUPSECNAMES = sfsem sfsstatistics figures tables
SUPSECS := $(addprefix $(SUPDIR)/, $(addsuffix .tex, $(SUPSECNAMES)))

SUPFIGNAMES = human_depth.png impala_depth.png yri_loglik.png masaimara_loglik.png f2.png stop.png stream.png bench_masaimara.png sim_truth.png sim_realsfs.png sim_winsfs.png sim_peak_truth.png sim_peak_realsfs.png sim_peak_winsfs.png
SUPFIGS := $(addprefix $(FIGDIR)/, $(SUPFIGNAMES))

SUPTABS := $(TABDIR)/data_test.tex

.PHONY: all clean view

all: $(MAIN) $(SUPPLEMENTARY) $(BOTH)

clean:
	rm -r $(TMPDIR) && rm -r $(OUTDIR)

view:
	xdg-open $(BOTH)

$(MAIN): $(DEPS) $(SECS) $(FIGS) $(TABS) $(SUPPLEMENTARY)

$(SUPPLEMENTARY): $(DEPS) $(SUPFIGS) $(SUPSECS) $(SUPTABS)

$(BOTH): $(MAIN) $(SUPPLEMENTARY)

$(OUTDIR)/%.pdf: %.tex
	latexmk -pdf -emulate-aux-dir -auxdir=$(TMPDIR) -outdir=$(OUTDIR) $<

$(FIGDIR)/%.pdf: $(TIKZDIR)/%.tex
	latexmk -pdf -emulate-aux-dir -auxdir=$(TMPDIR) -outdir=$(FIGDIR) $<

$(addprefix $(FIGDIR)/sim_, $(addsuffix .png, truth realsfs winsfs loglik)) &: $(RDIR)/sim.R $(FIGDEPS)
	Rscript $<

$(addprefix $(FIGDIR)/peak_, $(addsuffix .png, truth realsfs winsfs)) &: $(RDIR)/peak.R $(FIGDEPS)
	Rscript $<

$(FIGDIR)/%.png: $(RDIR)/%.R $(FIGDEPS)
	Rscript $<

$(TABDIR)/data_train.tex $(TABDIR)/data_test.txt: $(RDIR)/tables.R
	Rscript $<
