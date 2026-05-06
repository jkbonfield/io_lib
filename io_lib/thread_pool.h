/*
 * Copyright (c) 2013 Genome Research Ltd.
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
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 *
 * This file implements a thread pool for multi-threading applications.
 * It consists of two distinct interfaces: thread pools an results queues.
 *
 * The pool of threads is given a function pointer and void* data to pass in.
 * This means the pool can run jobs of multiple types, albeit first come
 * first served with no job scheduling.
 *
 * Upon completion, the return value from the function pointer is added to
 * a results queue. We may have multiple queues in use for the one pool.
 *
 * An example: reading from BAM and writing to CRAM with 10 threads. We'll
 * have a pool of 10 threads and two results queues holding decoded BAM blocks
 * and encoded CRAM blocks respectively.
 */

#ifndef _THREAD_POOL_H_
#define _THREAD_POOL_H_

// Convert hts thread pool to io_lib nomenclature.
// This is an easy one to one mapping as htslib's implementation started
// off here in io_lib.
#include <htslib/thread_pool.h>
typedef hts_tpool t_pool;
typedef hts_tpool_process t_results_queue;
typedef hts_tpool_result t_pool_result;

#define t_pool_init(q,t)             hts_tpool_init((t))
#define t_pool_dispatch              hts_tpool_dispatch
#define t_pool_dispatch2             hts_tpool_dispatch2
#define t_pool_flush                 hts_tpool_process_flush
#define t_pool_next_result           hts_tpool_next_result
#define t_pool_next_result_wait      hts_tpool_next_result_wait
#define t_pool_delete_result         hts_tpool_delete_result
#define t_pool_destroy(p,k)          hts_tpool_destroy((p))
#define t_pool_results_queue_empty   hts_tpool_process_empty
#define t_pool_delete_result         hts_tpool_delete_result
#define t_results_queue_init         hts_tpool_process_init
#define t_results_queue_destroy      hts_tpool_process_destroy
#define t_pool_results_queue_empty   hts_tpool_process_empty
#define t_pool_results_queue_len     hts_tpool_process_len
#define t_pool_results_queue_sz      hts_tpool_process_sz

#endif /* _THREAD_POOL_H_ */
