/*
 * Copyright (c) 2021 German Tischler-Hoehle
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

/**
 * test program reading an alignment file from stdin via scram_open/scram_get_seq
 * interface and writing it as a CRAM file on stdout
 *
 * The header of the alignment file needs to have at least the SN, LN, M5 and UR
 * attribute in each @SQ line
 *
 * Each reference sequence must be availabe in the reference cache given via the
 * REF_CACHE environment variable
 *
 * Example call:
 *
 * REF_CACHE=$HOME/refcache cram_bambam_interface <in.bam >out.cram
 **/

#include <io_lib/os.h>

#if defined(CRAM_IO_CUSTOM_BUFFERING)

#include <io_lib/cram_bambam.h>
#include <io_lib/scram.h>
#include <assert.h>
#include <string.h>
#include <pthread.h>

/**
 * data structure containing a semaphore
 **/
typedef struct _Semaphore
{
    pthread_mutex_t mutex;
    int mutexinit;
    pthread_cond_t cond;
    int condinit;
    size_t volatile c;
} Semaphore;

/**
 * deallocate a semaphore, returns NULL
 **/
Semaphore * deallocateSemaphore(Semaphore * sem)
{
    if ( sem )
    {
        if ( sem->mutexinit )
            pthread_mutex_destroy(&(sem->mutex));
        if ( sem->condinit )
            pthread_cond_destroy(&(sem->cond));

        free(sem);
    }

    return NULL;
}

/**
 * allocates a semaphore, returns NULL on failure
 **/
Semaphore * allocateSemaphore()
{
    Semaphore * sem = NULL;

    sem = (Semaphore *)malloc(sizeof(Semaphore));

    if ( ! sem )
        return deallocateSemaphore(sem);

    sem->mutexinit = 0;
    sem->condinit = 0;
    sem->c = 0;

    if ( pthread_mutex_init(&(sem->mutex),NULL) != 0 )
        return deallocateSemaphore(sem);

    sem->mutexinit = 1;

    if ( pthread_cond_init(&(sem->cond),NULL) != 0 )
        return deallocateSemaphore(sem);

    sem->condinit = 1;

    return sem;
}

/**
 * post operation for semaphore, returns 0 on success, nonzero otherwise
 **/
int semaphorePost(Semaphore * sem)
{
    if ( pthread_mutex_lock(&(sem->mutex)) )
        return -1;

    sem->c += 1;

    if ( pthread_cond_broadcast(&(sem->cond)) )
    {
        pthread_mutex_unlock(&(sem->mutex));
        return -1;
    }

    pthread_mutex_unlock(&(sem->mutex));

    return 0;
}

/**
 * wait operation for semaphore, returns 0 on success, nonzero otherwise
 **/
int semaphoreWait(Semaphore *sem)
{
    if ( pthread_mutex_lock(&(sem->mutex)) != 0 )
    {
        fprintf(stderr,"Failed to lock semaphore mutex\n");
        return -1;
    }

    while ( 1 )
    {
        if ( sem->c )
        {
            sem->c -= 1;
            pthread_mutex_unlock(&(sem->mutex));
            return 1;
        }

        if ( pthread_cond_wait(&(sem->cond),&(sem->mutex)) != 0 )
        {
            fprintf(stderr,"Failed to cond wait\n");
            return -1;
        }
    }
}

/**
 * data block
 **/
typedef struct _DataBlock
{
    /* memory block */
    char * cmem;
    /* size of memory block */
    size_t cmem_n;
    /* number of elements used in memory block */
    size_t cmem_f;

    /* start pointers to records in cmem */
    char ** pmem;
    /* size of pmem */
    size_t pmem_n;
    /* number of elements used in pmem */
    size_t pmem_f;

    /* array of element sizes stored by cmem */
    size_t * blocklen;
    size_t blocklen_n;
    size_t blocklen_f;

    /* array of element offsets into cmem */
    size_t * blockoff;
    size_t blockoff_n;
    size_t blockoff_f;
} DataBlock;

/* adjust pointers in pmem */
static void adjustDataBlock(DataBlock * DB)
{
    for ( size_t i = 0; i < DB->blocklen_f; ++i )
    {
        // fprintf(stderr,"offset %lu\n", DB->blockoff[i]);
        DB->pmem[i] = DB->cmem + DB->blockoff[i];
    }
}

/**
 * allocate a data block, returns NULL on failure
 **/
static DataBlock * allocateDataBlock()
{
    DataBlock * DB = NULL;

    DB = (DataBlock *)malloc(sizeof(DataBlock));

    if ( ! DB )
        return NULL;

    DB->cmem = NULL;
    DB->cmem_n = 0;
    DB->cmem_f = 0;

    DB->pmem = NULL;
    DB->pmem_n = 0;
    DB->pmem_f = 0;

    DB->blocklen = NULL;
    DB->blocklen_n = 0;
    DB->blocklen_f = 0;

    DB->blockoff = NULL;
    DB->blockoff_n = 0;
    DB->blockoff_f = 0;

    return DB;
}

/**
 * free a data block
 **/
static void freeDataBlock(DataBlock * DB)
{
    if ( DB )
    {
        free(DB->cmem);
        DB->cmem = NULL;
        free(DB->pmem);
        DB->pmem = NULL;
        free(DB->blocklen);
        DB->blocklen = NULL;
        free(DB->blockoff);
        DB->blockoff = NULL;
    }

    free(DB);
}

/**
 * increase size of cmem block by at least one
 **/
static int dataBlockBumpCMEM(DataBlock * DB)
{
    size_t const nn = DB->cmem_n ? (2*DB->cmem_n) : 1;

    char * nmem = malloc(nn);

    if ( ! nmem )
        return -1;

    memcpy(nmem,DB->cmem,DB->cmem_f);
    free(DB->cmem);
    DB->cmem = nmem;
    DB->cmem_n = nn;

    return 0;
}

/**
 * bump size of a memory block
 *
 * @param mem pointer to memory block pointer
 * @param pointer to size of block
 * @param pointer to number of elements used in block
 **/
static int bumpSize(size_t ** mem, size_t * n, size_t * f)
{
    size_t const nn = *n ? (2*(*n)) : 1;

    size_t * nmem = (size_t *)malloc(nn*sizeof(size_t));

    if ( ! nmem )
        return -1;

    memcpy(nmem,*mem,*f * sizeof(size_t));
    free(*mem);

    *mem = nmem;
    *n = nn;

    return 0;
}

/**
 * same as bumpSize but for array of pointers
 **/
static int bumpPtr(char *** mem, size_t * n, size_t * f)
{
    size_t const nn = *n ? (2*(*n)) : 1;

    char ** nmem = (char **)malloc(nn*sizeof(char *));

    if ( ! nmem )
        return -1;

    memcpy(nmem,*mem,*f * sizeof(char *));
    free(*mem);

    *mem = nmem;
    *n = nn;

    return 0;
}

/**
 * push alignment record into data block
 **/
static int dataBlockPush(DataBlock * DB, char const * const data, size_t const n)
{
    size_t const offset = DB->cmem_f;
    size_t const numsize = sizeof(uint32_t);
    size_t i;

    /* make sure we have sufficient space */
    while ( DB->cmem_f + n + numsize > DB->cmem_n )
        if ( dataBlockBumpCMEM(DB) )
            return -1;
    while ( DB->pmem_f + 1 > DB->pmem_n )
        if ( bumpPtr(&(DB->pmem),&(DB->pmem_n),&(DB->pmem_f)) )
            return -1;
    while ( DB->blocklen_f + 1 > DB->blocklen_n )
        if ( bumpSize(&(DB->blocklen),&(DB->blocklen_n),&(DB->blocklen_f)) )
            return -1;
    while ( DB->blockoff_f + 1 > DB->blockoff_n )
        if ( bumpSize(&(DB->blockoff),&(DB->blockoff_n),&(DB->blockoff_f)) )
            return -1;

    assert ( DB->cmem_f + n + numsize <= DB->cmem_n );
    assert ( DB->blocklen_f + 1 <= DB->blocklen_n );
    assert ( DB->blockoff_f + 1 <= DB->blockoff_n );

    /* store size of record as little endian number */
    for ( i = 0; i < numsize; ++i )
    {
        int const shift = i*CHAR_BIT;
        DB->cmem[DB->cmem_f++] = (n >> shift) & 0xFF;
    }

    /* store actual data */
    memcpy(DB->cmem + DB->cmem_f, data, n);

    /* update block meta data */
    DB->cmem_f += n;
    DB->pmem[DB->pmem_f++]         = NULL;
    DB->blocklen[DB->blocklen_f++] = n;
    DB->blockoff[DB->blockoff_f++] = offset;

    return 0;
}

/**
 * number pair
 **/
typedef struct _InBlockIdPair
{
    size_t inblockid;
    size_t blockid;
} InBlockIdPair;

/**
 * encoding control data used as userdata in encoding
 **/
typedef struct _ControlData
{
    /* semaphore for handling data blocks */
    Semaphore * dataBlockSemaphore;
    /* pointers to data blocks */
    DataBlock ** inputBlocks;
    /* number of data blocks allocated */
    size_t numinputblocks;

    /* input block free list (storing indexes into inputBlocks) */
    size_t volatile * inputBlocksFreeList;
    /* number of elements currently in inputBlockFreeList */
    size_t volatile inputBlocksFreeList_f;

    /* (inblockid,datablockid) pairs for blocks currently active in encoding */
    InBlockIdPair volatile * inBlockIdPair;
    /* number of entries currently held in inBlockIdPair */
    size_t volatile inBlockIdPair_f;

    /* lock for lists */
    pthread_mutex_t inputBlocksFreeList_f_lock;
    /* nonzero if mutex has been initialized */
    int inputBlocksFreeList_f_lock_init;

    /* next output block (not currently used) */
    size_t nextoutblock;

    /* output file */
    FILE * outfile;
} ControlData;

#define RUN_CLEANUP_IF_FAILURE(flag) { if ( flag ) { rv = EXIT_FAILURE; goto cleanup; } }

/**
 * free control data structure
 **/
ControlData * freeControlData(ControlData * CD)
{
    if ( CD )
    {
        if ( CD->dataBlockSemaphore )
        {
            deallocateSemaphore(CD->dataBlockSemaphore);
            CD->dataBlockSemaphore = NULL;
        }
        if ( CD->inputBlocks )
        {
            for ( size_t i = 0; i < CD->numinputblocks; ++i )
                if ( CD->inputBlocks[i] )
                {
                    freeDataBlock(CD->inputBlocks[i]);
                    CD->inputBlocks[i] = NULL;
                }

            free(CD->inputBlocks);
        }
        if ( CD->inputBlocksFreeList_f_lock_init )
        {
            pthread_mutex_destroy(&(CD->inputBlocksFreeList_f_lock));
        }
        if ( CD->inputBlocksFreeList )
        {
            free((size_t *)(CD->inputBlocksFreeList));
            CD->inputBlocksFreeList = NULL;
        }
        if ( CD->inBlockIdPair )
        {
            free((InBlockIdPair *)(CD->inBlockIdPair));
            CD->inBlockIdPair = NULL;
        }
        free(CD);
    }

    return NULL;
}

/**
 * allocate control data structure
 **/
ControlData * allocateControlData(size_t const rnuminputblocks, FILE * outfile)
{
    ControlData * CD = NULL;
    size_t i;

    CD = (ControlData*)malloc(sizeof(ControlData));

    if ( ! CD )
        return freeControlData(CD);

    CD->dataBlockSemaphore = NULL;
    CD->inputBlocks = NULL;
    CD->nextoutblock = 0;
    CD->numinputblocks = 0;
    CD->inBlockIdPair = 0;

    CD->inputBlocksFreeList_f_lock_init = 0;
    CD->inputBlocksFreeList_f = 0;
    CD->inputBlocksFreeList = NULL;
    CD->inBlockIdPair_f =0;

    CD->outfile = outfile;

    CD->dataBlockSemaphore = allocateSemaphore();

    if ( ! CD->dataBlockSemaphore )
        return freeControlData(CD);

    CD->inputBlocks = (DataBlock **)malloc(rnuminputblocks * sizeof(DataBlock *));
    if ( ! CD->inputBlocks )
        return freeControlData(CD);

    CD->numinputblocks = rnuminputblocks;
    for ( size_t i = 0; i < CD->numinputblocks; ++i )
        CD->inputBlocks[i] = NULL;
    for ( size_t i = 0; i < CD->numinputblocks; ++i )
        if ( (CD->inputBlocks[i] = allocateDataBlock()) == NULL )
            return freeControlData(CD);

    for ( size_t i = 0; i < CD->numinputblocks; ++i )
        if ( semaphorePost(CD->dataBlockSemaphore) != 0 )
            return freeControlData(CD);

    CD->inputBlocksFreeList = (size_t *)malloc(rnuminputblocks * sizeof(size_t));

    if ( ! CD->inputBlocksFreeList )
        return freeControlData(CD);

    for ( i = 0; i < CD->numinputblocks; ++i )
        CD->inputBlocksFreeList[CD->inputBlocksFreeList_f++] = i;

    CD->inBlockIdPair = (InBlockIdPair *)malloc(rnuminputblocks * sizeof(InBlockIdPair));

    if ( ! CD->inBlockIdPair )
        return freeControlData(CD);

    if ( pthread_mutex_init(&(CD->inputBlocksFreeList_f_lock),NULL) != 0 )
        return freeControlData(CD);

    return CD;
}

DataBlock * getDataBlock(ControlData * CD, size_t const inblockid)
{
    size_t i;

    /* wait until there is a free block */
    if ( semaphoreWait(CD->dataBlockSemaphore) != 1 )
    {
        fprintf(stderr,"getDataBlock: failed semaphore wait\n");
        pthread_mutex_unlock(&(CD->inputBlocksFreeList_f_lock));
        return NULL;
    }

    /* get lock for free lists */
    if ( pthread_mutex_lock(&(CD->inputBlocksFreeList_f_lock)) != 0 )
    {
        fprintf(stderr,"getDataBlock: failed to lock mutex\n");
        return NULL;
    }

    /* list should not be empty */
    assert ( CD->inputBlocksFreeList_f );
    i = CD->inputBlocksFreeList[--CD->inputBlocksFreeList_f];

    /* store (inblockid,blockid) pair */
    InBlockIdPair IBIP;
    IBIP.inblockid = inblockid;
    IBIP.blockid = i;
    CD->inBlockIdPair[CD->inBlockIdPair_f++] = IBIP;

    /* release lock for free lists */
    pthread_mutex_unlock(&(CD->inputBlocksFreeList_f_lock));

    /* get data block pointer */
    DataBlock * block = CD->inputBlocks[i];

    /* reset data block */
    block->cmem_f = 0;
    block->pmem_f = 0;
    block->blocklen_f = 0;
    block->blockoff_f = 0;

    return block;
}

/**
 * return data block specified by inblockid
 **/
void putDataBlock(ControlData * CD, size_t const inblockid)
{
    /* get lock */
    if ( pthread_mutex_lock(&(CD->inputBlocksFreeList_f_lock)) != 0 )
        return;

    size_t blockid = 0;
    size_t i = 0;
    size_t o = 0;
    int found = 0;

    /* look for (inblock,blockid) pair in list and erase it */
    for ( i = 0; i < CD->inBlockIdPair_f; ++i )
        if ( CD->inBlockIdPair[i].inblockid == inblockid )
        {
            blockid = CD->inBlockIdPair[i].blockid;
            found = 1;
        }
        else
        {
            CD->inBlockIdPair[o++] = CD->inBlockIdPair[i];
        }

    /* set new list size */
    CD->inBlockIdPair_f = o;

    assert ( found );

    CD->inputBlocksFreeList[CD->inputBlocksFreeList_f++] = blockid;

    pthread_mutex_unlock(&(CD->inputBlocksFreeList_f_lock));

    semaphorePost(CD->dataBlockSemaphore);
}

void enque_compression_work_package_function(void *userdata, void *workpackage)
{
    int const r = cram_process_work_package(workpackage);
    assert ( r == 0 );
}

void compression_work_package_finished(void *userdata, size_t const inblockid, int const final)
{
    ControlData * CD = (ControlData*)userdata;
    putDataBlock(CD,inblockid);
}

void write_function(void *userdata, ssize_t const inblockid, size_t const outblockid, char const *data, size_t const n, cram_data_write_block_type const blocktype)
{
    ControlData * CD = (ControlData*)userdata;
    /* fprintf(stderr,"write %lu\n", (unsigned long) n); */
    ssize_t const w = fwrite(data,n,1,CD->outfile);
    assert ( w == 1 );
}

int main(int argc, char *argv[])
{
    /* maximum input block size in alignment records */
    static unsigned int maxblocksize = 16*1024;
    /* number of input blocks */
    static unsigned int numinputblocks = 4;
    scram_fd * infd = NULL;
    ControlData * CD = NULL;
    int rv = EXIT_SUCCESS;
    SAM_hdr * hdr = NULL;
    SAM_hdr * hdrdup = NULL;
    int hdrlen = -1;
    char const * hdrtext = NULL;
    void * cram_encoder = NULL;
    bam_seq_t *bsp = NULL;
    size_t encoded = 0;
    size_t lastprint = 0;
    int const printshift = 20;
    char const * infn = NULL;
    char const * outfn = NULL;
    FILE * outfile = NULL;

    if ( !(2 < argc) )
    {
        fprintf(stderr,"usage %s <in> <out>\n", argv[0]);
        return EXIT_FAILURE;
    }

    infn = argv[1];
    outfn = argv[2];

    infd = scram_open(infn,"r");

    RUN_CLEANUP_IF_FAILURE(!infd);

    outfile = fopen(outfn,"w");

    RUN_CLEANUP_IF_FAILURE(!outfile);

    CD = allocateControlData(numinputblocks,outfile);

    RUN_CLEANUP_IF_FAILURE(!CD);

    hdr = scram_get_header(infd);

    RUN_CLEANUP_IF_FAILURE(!hdr);

    hdrdup = sam_hdr_dup(hdr);

    RUN_CLEANUP_IF_FAILURE(!hdrdup);

    hdrlen = sam_hdr_length(hdrdup);

    RUN_CLEANUP_IF_FAILURE(hdrlen < 0);

    hdrtext = sam_hdr_str(hdrdup);

    RUN_CLEANUP_IF_FAILURE(!hdrtext);

    cram_encoder = cram_allocate_encoder(CD,hdrtext,hdrlen,write_function);

    RUN_CLEANUP_IF_FAILURE(!cram_encoder);

    for ( size_t inblockid = 0; !scram_eof(infd); ++inblockid )
    {
        size_t i = 0;
        DataBlock * DB = getDataBlock(CD,inblockid);

        assert ( DB );

        /* read data (up to maxblocksize alignment records) */
        while ( i < maxblocksize && scram_get_seq(infd,&bsp) == 0 )
        {
            char const * const data = (char const *)(&(bsp->ref));
            size_t const len = bam_blk_size(bsp);
            int const r = dataBlockPush(DB, data, len);
            assert ( r == 0 );
            i += 1;
        }

        adjustDataBlock(DB);

	char const *block[]        = { DB->cmem };
	size_t const blocksize     = DB->cmem_f;
	size_t const blockelements = DB->blocklen_f;
	
        int const rq = cram_enque_compression_block(
            CD,
            cram_encoder,
            inblockid,
            block,
            &blocksize,
            &blockelements,
            1,
            scram_eof(infd)?1:0, /* final */
            enque_compression_work_package_function,
            write_function,
            compression_work_package_finished
        );

        assert ( rq == 0 );

        encoded += i;

        if ( (encoded >> printshift) != (lastprint >> printshift) )
        {
            fprintf(stderr,"%lu\n",(unsigned long)encoded);
            lastprint = encoded;
        }
    }

    fprintf(stderr,"%lu\n",(unsigned long)encoded);

    cleanup:
    if ( cram_encoder )
        cram_deallocate_encoder(cram_encoder);
    if ( hdrdup )
        sam_hdr_free(hdrdup);
    if ( infd )
        scram_close(infd);
    if ( outfile )
        fclose(outfile);
    freeControlData(CD);

    return rv;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif
