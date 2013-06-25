## Makefile.am -- Process this file with automake to produce Makefile.in

EXTRA_DIST              = $(TESTS) data compare_sam.pl generate_data.pl
MAINTAINERCLEANFILES    = Makefile.in

test_outdir              = test.out

TESTS_ENVIRONMENT       = \
			top_builddir=$(top_builddir) \
			outdir=$(test_outdir)

TESTS                   = \
                        srf2fasta.test \
                        srf2fastq.test \
			srf_info.test \
			srf_filter.test \
			srf_index.test \
			scram.test \
			scram_mt.test

dist-hook:
	rm -rf `find $(distdir)/data -name .svn`
	rm -f $(distdir)/data/ce#sorted.bam $(distdir)/data/ce#unsorted.bam
	rm -f $(distdir)/data/*~

distclean-local:
	-rm -rf $(test_outdir)

