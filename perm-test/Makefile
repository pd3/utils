# Copyright (C) 2018-2020 Genome Research Ltd.
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


# Create the config.mk file manually to override the default settings
ifneq "$(wildcard config.mk)" ""
include config.mk
endif

HTSDIR ?= ../htslib
include $(HTSDIR)/htslib_static.mk

CC     ?= gcc
CFLAGS ?= -Wall -g -O2 -I ../libs -I $(HTSDIR)
LDLIBS ?= ../libs/regidx.o ../libs/dist.o $(HTSDIR)/libhts.a $(HTSLIB_static_LIBS)
BIN     = perm-test recurrence-test


all: $(BIN)

perm-test.o: ../libs/version.h ../libs/dist.h perm-test.c

recurrence-test.o: ../libs/version.h recurrence-test.c

perm-test: perm-test.o
	$(CC) $(LDFLAGS) perm-test.o $(LDLIBS) -o $@

recurrence-test: recurrence-test.o
	$(CC) $(LDFLAGS) recurrence-test.o $(LDLIBS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $<

test: $(BIN)
	./test/test.pl

clean:
	rm -f $(BIN) *.o *.a *~

force:

.PHONY: all clean test

