bin_PROGRAMS = 
bin_SCRIPTS =
man_MANS = 
doc_DATA = bammds_flowchart.pdf Gus.Subsample.bam HDGP.Reference.Subsample.txt TUTORIAL

all: bammds_flowchart.pdf

bammds_flowchart.pdf:  bammds_flowchart.dot
	dot -Tpdf -o bammds_flowchart.pdf bammds_flowchart.dot

DISTCLEANFILES = bammds_flowchart.pdf
EXTRA_DIST = bammds_flowchart.dot \
	bammds_flowchart.pdf \
	Gus.Subsample.bam \
	HDGP.Reference.Subsample.txt \
	TUTORIAL
