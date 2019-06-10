# Copyright (C) 2019 Genome Research Ltd.
#
# Author: Petr Danecek <pd3@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


# Edit the location of HTSlib here or run as HTSDIR=/path/to/htslib make
export HTSDIR = ../htslib


TARGETS = annot-regs

all:
	@$(MAKE) libs/version.h; \
    ln -fs $(HTSDIR) htslib; \
    $(MAKE) -C libs; \
	for T in $(TARGETS) ; do $(MAKE) -C $$T; done

clean:
	@$(MAKE) -C libs clean; \
    for T in $(TARGETS) ; do $(MAKE) -C $$T clean; done

test:
	@for T in $(TARGETS) ; do $(MAKE) -C $$T test; done


PACKAGE_VERSION := 1.0.$(shell git describe --always --dirty)
libs/version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))

libs/version.h:
	echo '#define ANNOT_REGS_VERSION "$(PACKAGE_VERSION)"' > $@




