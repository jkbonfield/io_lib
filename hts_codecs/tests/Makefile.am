# Copyright (c) 2019 Genome Research Ltd.
# Author(s): James Bonfield
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# 
#    1. Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
# 
#    2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
# 
#    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
#    Institute nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific
#    prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
# LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 

# Standalone test programs
bin_PROGRAMS = rans4x16pr tokenise_name3 arith_dynamic rans4x8 rans4x16pr fqzcomp_qual

LDADD = $(top_builddir)/hts_codecs/libhts-codecs.la -lm
CPPFLAGS = -I$(top_srcdir)

fqzcomp_qual_SOURCES = fqzcomp_qual_test.c
rans4x8_SOURCES = rANS_static_test.c
rans4x16pr_SOURCES = rANS_static4x16pr_test.c
arith_dynamic_SOURCES = arith_dynamic_test.c
tokenise_name3_SOURCES = tokenise_name3_test.c

TESTS = \
	rans4x8.test \
	rans4x16.test \
	arith.test \
	tok3.test \
	fqzcomp.test

EXTRA_DIST = $(TESTS) dat names

test_outdir = test.out

distclean-local:
	-rm -rf $(test_outdir)
