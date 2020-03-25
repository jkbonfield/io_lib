/*
 * Copyright (c) 2015 Genome Research Ltd.
 * Author(s): James Bonfield
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 * 
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *    Institute nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific
 *    prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


/*
 * This code contains an interface for Biobambam and Libmaus to use
 * when encoding CRAM files.
 *
 * Libmaus has its own thread-pool system.  We construct work packages
 * to be dispatched to libmaus.
 *
 * The basic model is that we have a single encoder context, created /
 * destroyed by cram_allocate_encoder() and cram_deallocate_encoder().
 * The context can just be a cram_fd pointer.
 *
 * This context is then passed along with a block of uncompressed BAM
 * records to cram_enque_compression_block(), running in the top thread.
 * This procedure creates a work package and puts it on the work queue
 * passed in. 
 *
 * The cram_process_work_package() function is then called per work
 * package, within 1 or more threads.  Within this function it will
 * call the libmaus supplied write function and when finished will
 * call the libmaus supplied work-finished function.  The write
 * function takes the place of the CRAM_IO_PUTC and CRAM_IO_WRITE
 * functions already defined in cram_io.c.
 */
#if ! defined(IO_LIB_CRAM_BAMBAM_H)
#define IO_LIB_CRAM_BAMBAM_H

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "io_lib/cram.h"

//-----------------------------------------------------------------------------
// Public interface

typedef enum cram_data_write_block_type {
    cram_data_write_block_type_internal,
    cram_data_write_block_type_block_final,
    cram_data_write_block_type_file_final
} cram_data_write_block_type;

// Enqueue a package of work to compress a CRAM slice.
typedef void
(*cram_enque_compression_work_package_function_t)(void *userdata,
						  void *workpackage);

// Callback to indicate the block has been compressed
typedef void
(*cram_compression_work_package_finished_t)(void *userdata,
					    size_t const inblockid,
					    int const final);

// Write function for compressed blocks, provided by libmaus.
// Inblockid is the same as supplied by the dispatcher.
// Outblockid increments from 0 per unique inblockid.
typedef void
(*cram_data_write_function_t)(void *userdata,
			      ssize_t const inblockid,
			      size_t const outblockid,
			      char const *data,
			      size_t const n,
			      cram_data_write_block_type const blocktype);


// Temporary copy from biobambam (BSD licence verbally accepted) to
// help validate the interface via the compiler.
extern void *cram_allocate_encoder(void *userdata,
				   char const *sam_header,
				   size_t const sam_headerlength,
				   cram_data_write_function_t write_func);
extern void cram_deallocate_encoder(void *context);
extern int cram_enque_compression_block(
	void *userdata,
	void *context,
	size_t const inblockid,
	char const **block,
	size_t const *blocksize,
	size_t const *blockelements,
	size_t const numblocks,
	int const final,
	cram_enque_compression_work_package_function_t workenqueuefunction,
	cram_data_write_function_t writefunction,
	cram_compression_work_package_finished_t workfinishedfunction);
extern int cram_process_work_package(void *workpackage);
extern cram_fd * cram_encoder_get_fd(void *context);
#endif
