# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXBUILD   = sphinx-build
VERSION = $(shell cat ../version.txt | grep -o -E [0-9]+.[0-9]+)
HTMLBUILDDIR  = ./build/html

# User-friendly check for sphinx-build
ifeq ($(shell which $(SPHINXBUILD) >/dev/null 2>&1; echo $$?), 1)
$(error The '$(SPHINXBUILD)' command was not found. Make sure you have Sphinx installed, then set the SPHINXBUILD environment variable to point to the full path of the '$(SPHINXBUILD)' executable. Alternatively you can add the directory with the executable to your PATH. If you don't have Sphinx installed, grab it from http://sphinx-doc.org/)
endif

# Internal variables.
ALLSPHINXOPTS   = -n -d $(HTMLBUILDDIR)/doctrees source

.PHONY: all help clean html

all: html

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"

clean:
	rm -rf $(HTMLBUILDDIR)

html:
	mkdir -p $(HTMLBUILDDIR)
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(HTMLBUILDDIR)
	@echo
	@echo "Build finished. The HTML pages are in $(HTMLBUILDDIR)."

view:
	xdg-open $(HTMLBUILDDIR)/index.html
