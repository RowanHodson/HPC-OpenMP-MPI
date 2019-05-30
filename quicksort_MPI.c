/*
 * Parallel Quicksort implementation using MPI
 *
 * Author: Krzysztof Voss [shobbo@gmail.com]
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

int NSIZE;

typedef struct {
    int numtasks;
    int rank;

    int pivot;
    int *set;
    int *buf;
    int size;
    int lim;

    int gl, gh;
} mdata;

#define TEST
    



/////////////////////////////////
static int
debug_tables(mdata *md) {
    int p;

    printf("%d:\t\t\tdebug_tables()\n", md->rank);

    //  printf("%d:\t\t\t\tset: ", md->rank);
    //  for ( p = 0; p < md->size; p++)
    //      printf("%d ", md->set[p]);
    //  printf("\n");
    fflush(stdout);

    printf("%d:\t\t\t\tpivot: %d: , size: %d, lim: %d, gl: %d, gh: %d\n", md->rank, md->pivot, md->size, md->lim, md->gl, md->gh);
    fflush(stdout);

    return 0;
}



/////////////////////////////////
static int
free_memory(mdata* md) {


    free(md->set);
    md->set = NULL;

    return 0;
}

/////////////////////////////////
static int
gather_values(mdata *md) {
    int i;

    if (md->rank == 0) {
        int osize, offset;
        MPI_Status status;

        md->buf = (int*) malloc( NSIZE * sizeof(int) );
        /* because the lower the rank is the higher number we have */

        if ( md->size > 0 ) {
            memcpy( &md->buf[NSIZE - md->size], md->set, md->size * sizeof(int));
            offset = NSIZE - md->size;
        } else offset = NSIZE;

        free(md->set);
        md->set = md->buf;
        md->size = NSIZE;


        for ( i = 1; i < md->numtasks; i++) {
            MPI_Recv(&osize, 1, MPI_INT,  i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            offset -= osize;
            if ( offset < 0 ) {

                return 1;
            }
            MPI_Recv(md->set + offset, osize, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        }
   
    } else {
        int tag = 0;

        MPI_Send(&md->size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        MPI_Send(md->set, md->size, MPI_INT, 0, tag, MPI_COMM_WORLD);

    }

    return 0;
}

/////////////////////////////////
static int
load_values(mdata *md) {
    FILE *fp;
    int size;

    if ( md->rank != 0 ) return 0;

    fp = fopen("list-to-sort.txt", "rb");
    if ( fp == NULL ) return 1;

    fseek(fp, 0, SEEK_END);
    size = ftell(fp);

    md->set = (int*)malloc(size*sizeof(int));

    fclose(fp);

    return 0;
}

/////////////////////////////////
/*
 * funkcja moze byc wywolywana przez wszystkie procesy
 * tylko odpowiednie ranki ja wywolaja
 */
static int
send_buf(mdata *md, int *buf, int size, int dest, int source) {
    int p;

    /* should we call this procedure? */
    /* we are either the source or the destination */
    if ( md->rank == source ) {
        int tag = 0;

        /* we need to send an information on how big buffer we are going to send */
        MPI_Send(&size, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);

        MPI_Send(buf, size, MPI_INT, dest, tag, MPI_COMM_WORLD);


        return 0;
    }
    else if ( md->rank == dest ) {
        MPI_Status status;

        /* receiving the buffer size and allocating memory */
        MPI_Recv(&md->size, 1, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        md->set = (int*)malloc(  md->size * sizeof(int));
        MPI_Recv(md->set, md->size, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        return 0;
    }

    return 0;
}

/////////////////////////////////
static int
distribute_values(mdata *md) {
    int i;
    int *buf;
    int size;
    int offset;

#ifdef DEBUG
    printf("%d:\tdistribute_values()\n", md->rank);
    fflush(stdout);
#endif
    size = floor( NSIZE / md->numtasks );
    offset = size + ( NSIZE - md->numtasks * size );
    if ( md->rank == 0 )
        md->size = offset;

    /* sending parts to other processes */
    for ( i = 1; i < md->numtasks; i++ ) {
        // send_buf(mdata *md, int *buf, int size, int dest, int source)
        send_buf(md, md->set+offset, size, i, 0);
        offset += size;
    }

    return 0;
}

/////////////////////////////////
static int
set_pivot(mdata *md) {
    int i;
    int msg;
    MPI_Status status;

    /* if we have the lowest rank in the group we are the leader */
    if ( md->rank == md->gl ) {
        /* first we have to check if we have any data to choose pivot from */
        if ( md->size > 0 ) {
            msg = md->pivot = md->set[0];
            /* sending others the pivot */
            for ( i = md->rank + 1; i < md->gh; i++ )
                MPI_Send(&msg, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            return 0;
        }
        else {
            /* now we ask others in the group if they have a pivot */
            for( i = md->rank + 1; i < md->gh; i++ ) {
                MPI_Send(&msg, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Recv(&msg, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                /* if we were sent back a pivot */
                if ( status.MPI_TAG == 1 ) {
                    md->pivot = msg;
                    /* ordering others to set the pivot */
                    for( i = md->rank + 1; i < md->gh; i++ )
                        MPI_Send(&msg, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
                    /* the job is done */
                    return 0;
                }
            }
            /* noone has it, so we inform all to return */
            for( i = md->rank + 1; i < md->gh; i++ )
                MPI_Send(&msg, 1, MPI_INT, i, 2, MPI_COMM_WORLD);

        }
    } else {
        while (1) {
            /* waiting for a pivot */
            MPI_Recv(&msg, 1, MPI_INT, md->gl, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if ( status.MPI_TAG == 1 ) {
                /* we have a pivot */
                md->pivot = msg;
                return 0;
            } else if ( status.MPI_TAG == 0 ) {
                /* the leader doesn't have any pivot */
                /* do we have any? */
                if ( md->size > 0 ) {
                    /* one has it */
                    msg = md->set[0];
                    MPI_Send(&msg, 1, MPI_INT, md->gl, 1, MPI_COMM_WORLD);
                } else {
                    /* one does not have it */
                    MPI_Send(&msg, 1, MPI_INT, md->gl, 0, MPI_COMM_WORLD);
                }
            } else {
                /* we don't have any pivots set or other value was sent us */
#ifdef DEBUG
                printf("%d: set_pivot[%d %d] -> Status: %d\n", md->rank, md->gl, md->gh, status.MPI_TAG);
                fflush(stdout);
#endif
                return 1;
            }
        }
    }

    return 1;
}

/////////////////////////////////
static int
divide_numbers(mdata *md) {

    int k, l;
    int *p, *q;
    int buf;

    /* if we have any data to process */
    if ( md->size > 0 ) {
        l = md->size;

        /* initializing the first q */
        if (md->set) q = &md->set[l-1];
        else return 1;  
    } else return 1;

    /* starting with the beggining of the table to be sorted */
    for ( k = 0; k < l; k++ ) {
        p = &md->set[k];
        if ( *p > md->pivot ) {
            /* starting with the end we check if every latter number is bigger */
            while ( k < l ) {
                q = &md->set[l-1];
                if ( *q > md->pivot ) {
                    l--;
                    continue;
                } else {
                    buf = *p;
                    *p = *q;
                    *q = buf;
                    break; 
                }
            }
        }
    }
    /* returning result to the data structure */
    md->lim = l;

    return 0;
}

/////////////////////////////////
static int
swap_numbers(mdata *md) {
    int i;
    int numswaps;


    /* number of swaps in my group */
    numswaps = (md->gh - md->gl)/2;
    /* for every pair */
    for ( i = 0; i < numswaps; i++ ) {

        /* only two processes are involved in exchanging their data */
        int p0, p1;
        MPI_Status status;
        int osize;
        int lsize;
        int tag;

        p0 = md->gl + i;
        p1 = md->gl + numswaps + i;
        if ( md->rank == p0 ){
            int xhs;
            tag = 0;
            xhs = md->size - md->lim;
            lsize = md->lim;

            /* sending the lower part of my set */
            /* first we send an information about the size of our part */
            MPI_Send(&lsize, 1, MPI_INT, p1, tag, MPI_COMM_WORLD);

            /* now we are receiving the same information */
            MPI_Recv(&osize, 1, MPI_INT, p1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            /* now we allocate all necessary memory */
            md->buf = (int*) malloc( (xhs + osize)*sizeof(int) );
            memcpy ( md->buf + osize, md->set + md->lim, xhs * sizeof(int));

            /* we send our data */
            MPI_Send(md->set, lsize, MPI_INT, p1, tag, MPI_COMM_WORLD);
            free(md->set);
            md->set = md->buf;
            md->size = xhs + osize;

            /* and receive new */
            MPI_Recv(md->set, osize, MPI_INT, p1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            md->lim = 0;      
        } else if ( md->rank == p1 ) {
            tag = 0;            
            lsize = md->size - md->lim;


            /* waiting for the info about the amount of incoming data */
            MPI_Recv(&osize, 1, MPI_INT, p0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            /* now we allocate all necessary memory */
            md->buf = (int*) malloc((md->lim + osize)*sizeof(int));
            memcpy(md->buf, md->set, md->lim * sizeof(int));

            /* now we send our corresponding information  and confirm having mem allocated */
            MPI_Send(&lsize, 1, MPI_INT, p0, tag, MPI_COMM_WORLD);

            /* now receiving actual data */

            MPI_Recv(md->buf + md->lim, osize, MPI_INT, p0, tag, MPI_COMM_WORLD, &status);

            /* and sending ours */
            MPI_Send(md->set + md->lim, lsize, MPI_INT, p0, tag, MPI_COMM_WORLD);

            free(md->set);
            md->set = md->buf;
            md->size = md->lim + osize;
            md->lim = md->size;

        }
    }

    return 0;
}

/////////////////////////////////
static int
comp(const void *p, const void *q) {

    return *(int *)p - *(int *)q;
}

/////////////////////////////////
static int
pqsort(mdata *md) {
    int tasks, tasks2;

    /* if my rank doesn't belong here */
    if ( md->gl > md->rank || md->rank >= md->gh )
        return 0;
    /* if we're not sharing this with anyone */
    if ( md->gl == (md->gh - 1) ) {

        qsort(md->set, md->size, sizeof(int), comp);
        return 0;
    }

    /* if we don't have any pivot. we don't have anything to sort */
    if ( set_pivot(md) ) return 1;

    /* dividing values locally to get the md->lim */
    divide_numbers(md);
    /* swapping numbers with the parter in group */
    swap_numbers(md);

    /* we need to save some values on stack */
    //tasks = (md->gh - md->gl);
    tasks2 = (md->gh - md->gl) / 2;

    /* sorting the upper half */
    md->gh = md->gl + tasks2;
    pqsort(md);

    /* sorting the lower half */
    md->gl = md->gh;
    md->gh = md->gh + tasks2;
    pqsort(md);

    return 0;
}

/////////////////////////////////
static int
gen_set(mdata *md) {
    int i;

    /* if we are the leader */
    if ( md->rank == 0 ) {
        md->set = (int*)malloc(NSIZE * sizeof(int));
        md->size = NSIZE;
        srand(time(NULL));
        /* initializing array with random numbers */
        for (i=0; i < NSIZE; i++)
            md->set[i] = 1 + (int) ( 1024.0 * 1024.0 * 1024.0 * rand() / (RAND_MAX + 1.0));
    }

    return 0;
}

/////////////////////////////////
/*
 * returns 1 if sorting went fine
 */
static int
sort_ok(mdata *md) {
    int s = md->size - 1;

    while ( s-- > 0 ) {
        if ( md->set[s] > md->set[s+1] )
            return 0;
    }
    return 1;
}

/////////////////////////////////
int
main(int argc, char **argv) {
    mdata md;
    if (argv[0] != NULL)
    {
        NSIZE = atoi(argv[1]);
    }
#ifdef TEST
    
    double etime;   /* elapsed time */
#endif
    md.size = 0;
    md.gl = 0;
    md.set = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &md.numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &md.rank);
    md.gh = md.numtasks;
    gen_set(&md);

#ifdef TEST
    MPI_Barrier(MPI_COMM_WORLD);
    etime = - MPI_Wtime();
#endif
    distribute_values(&md);
    pqsort(&md);
    gather_values(&md);
#ifdef TEST
    etime += MPI_Wtime();
    if ( md.rank == 0 ) {
        if ( sort_ok(&md) )
            printf("0: Elapsed time: %.3f [s], numproc: %d, datasize: %d, status: OK\n", etime, md.numtasks, NSIZE);
        else
            printf("0: Elapsed time: %.3f [s], numproc: %d, datasize: %d, status: FAILED\n", etime, md.numtasks, NSIZE);
    }
#endif
    free_memory(&md);

    MPI_Finalize();

    return 0;
}
