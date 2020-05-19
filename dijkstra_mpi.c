/* assert */
#include <assert.h>
/* INFINITY */
#include <math.h>
/* FILE, fopen, fclose, fscanf, rewind */
#include <stdio.h>
/* EXIT_SUCCESS, malloc, calloc, free */
#include <stdlib.h>
/* time, CLOCKS_PER_SEC */
#include <time.h>
#include <mpi.h>

#define ROWMJR(R,C,NR,NC) (R*NC+C)
#define COLMJR(R,C,NR,NC) (C*NR+R)
/* define access directions for matrices */
#define a(R,C) a[ROWMJR(R,C,ln,n)]
#define b(R,C) b[ROWMJR(R,C,nn,n)]
#define MAIN_PROCESS 0
#define SEND_NUM_TAG 0
#define SEND_OFFSET_TAG 1
#define SEND_ELEMNTS_TAG 2
#define SEND_WEIGHT_TAG 3
#define SEND_COUNTS_TAG 4
#define SEND_RESULT_TAG 5

static void calculate_offset(
    int ** offset, /* offset[i]: starting vertex that the ith compute node is responsible for */ 
    int ** nlocal, /* nlocal[i]: the number of vertices the ith compute node is responsible for */ 
    int npe, /* Number of compute nodes */ 
    int n /* Number of vertices in the graph */ 
    ) {
    int local, remainder;
    if (n < npe) {
        local = 1;
        *nlocal = malloc(npe * sizeof(int));
        int i = 0;
        for (; i < n; i++) {
            *(*nlocal + i) = local; /* nlocal[i] = local */ 
        }

        for (; i < npe; i++) {
            *(*nlocal + i) = 0; /* nlocal[i] = 0 */ 
        }
    } else {
        local = n / npe;
        remainder = n % npe;

        *nlocal = malloc(npe * sizeof(int));
        for (int i = 0; i < npe; i++) {
            *(*nlocal + i) = local;
        }
        *(*nlocal + npe - 1) += remainder;
    }

    *offset = malloc(npe * sizeof(int));
    for (int i = 0; i < npe; i++) {
        *(*offset + i) = 0; /* offset[i] = 0 */ 
        for (int j = 0; j < i; j++) {
            *(*offset + i) += *(*nlocal + j); /* offset[i] += nlocal[j] */ 
        }
    }
}

static void
read_file_and_send(
        char const * const filename,
        int * const np,
        float ** const ap,
        int npe,
        int ** offset,
        int ** nlocal
)
{
    int i, j, n, ret;
    FILE *fp = NULL;
    float *a = NULL;

    /* open the file */
    fp = fopen(filename, "r");
    assert(fp);

    /* get the number of nodes in the graph */
    ret = fscanf(fp, "%d", &n);
    assert(1 == ret);

    calculate_offset(offset, nlocal, npe, n);

    a = malloc(n * *(*nlocal) * sizeof(float));
    for (j = 0; j < *(* nlocal) * n; j++) {
        ret = fscanf(fp, "%f", &a[j]);
        assert(1 == ret);
    }
    *ap = a;

    for (i = 1; i < npe; i++) {
        a = malloc(n * *(*nlocal + i) * sizeof(float));
        MPI_Send(&n, 1, MPI_INTEGER, i, SEND_NUM_TAG, MPI_COMM_WORLD);
        MPI_Send(*offset, npe, MPI_INTEGER, i, SEND_OFFSET_TAG, MPI_COMM_WORLD);
        MPI_Send(*nlocal, npe, MPI_INTEGER, i, SEND_ELEMNTS_TAG, MPI_COMM_WORLD);

        for (j = 0; j < *(* nlocal + i) * n; j++) {
            ret = fscanf(fp, "%f", &a[j]);
            assert(1 == ret);
        }

        MPI_Send(&j, 1, MPI_INTEGER, i, SEND_COUNTS_TAG, MPI_COMM_WORLD);
        MPI_Send(a, j, MPI_FLOAT, i, SEND_WEIGHT_TAG, MPI_COMM_WORLD);
        free(a);
    }

    /* close file */
    ret = fclose(fp);
    assert(!ret);
    *np = n; 
    MPI_Barrier(MPI_COMM_WORLD); 
}

static void
recv_values_from_master(
        int * const np,
        float ** const ap,
        int npe,
        int ** offset,
        int ** nlocal
)
{
    int count, n;
    float *a = NULL; 
    MPI_Recv(&n, 1, MPI_INTEGER, 0, SEND_NUM_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    *offset = malloc(npe * sizeof(int));
    MPI_Recv(*offset, npe, MPI_INTEGER, MAIN_PROCESS, SEND_OFFSET_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    *nlocal = malloc(npe * sizeof(int)); 
    MPI_Recv(*nlocal, npe, MPI_INTEGER, MAIN_PROCESS, SEND_ELEMNTS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
    MPI_Recv(&count, 1, MPI_INTEGER, MAIN_PROCESS, SEND_COUNTS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    a = malloc(count * sizeof(float)); 
    MPI_Recv(a, count, MPI_FLOAT, MAIN_PROCESS, SEND_WEIGHT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
    *ap = a; 
    *np = n; 
    MPI_Barrier(MPI_COMM_WORLD); 
}

static void
dijkstra(
        int const s, // Source
        int const n, // Number of nodes
        float const * const a, // Adjacency matrix
        float **const lp, // Result
        int rank, 
        int * offset, /* offset[i]: starting vertex that the ith compute node is responsible for */ 
        int * nlocal, /* nlocal[i]: the number of vertices the ith compute node is responsible for */ 
        int npe /* Number of compute nodes */ 
)
{
    int i, j, source_node = 0; 
    struct float_int {
        float l;
        int u;
    } min;
    
    char * m; /* m[i] == 1: i is visited */
    float *l = NULL;
    float *local_result = NULL;

    m = calloc(n, sizeof(*m));
    assert(m);
    l = malloc(n * sizeof(float)); 
    assert(l); 
    local_result = malloc(n * sizeof(float)); 
    assert(local_result);  

    //TODO: Refactor to have source_node sent from main process to all other processes 
    //TODO: Refactor to have only specific offset and nlocal sent to all processes 
    for (i = 0; i < npe; i++) {
        /* Find the compute node that has the source vertex */ 
        if (s < offset[i]) {
            source_node = i - 1; 
            break; 
        }
    }

    /* If I (this compute node) has the source vertex */
    if (rank == source_node) {
        for (i = 0; i < n; i++) {
            /* Initialize distances to the distance from source to all other vertices */ 
            l[i] = a((s - offset[source_node]), i); 
        }
    }
    
    MPI_Bcast(l, n, MPI_FLOAT, source_node, MPI_COMM_WORLD); /* Broadcast from source node to all others */ 
    MPI_Barrier(MPI_COMM_WORLD); 

    m[s] = 1; /* source vertex visited */ 
    min.u = -1; 

    for (i = 1; i < n; i++) {
        min.l = INFINITY;
        for (j = 0; j < n; j++) {
            if (!m[j] && l[j] < min.l) {
                min.l = l[j];
                min.u = j;
            }
            local_result[j] = l[j]; 
        }

        m[min.u] = 1; 
        for (j = 0; j < nlocal[rank]; j++) {
            if (m[j + offset[rank]]) {
                /* If already visited */ 
                continue;
            }
            /* If a shorter path is found */ 
            if (a(j, min.u) + min.l < local_result[j + offset[rank]]) {
                local_result[j + offset[rank]] = a(j, min.u) + min.l;
            }
        }
        /* Send out local_result for MIN reduction. Final result is received in l */ 
        MPI_Allreduce(local_result, l, n, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD); 
        MPI_Barrier(MPI_COMM_WORLD); 
    }
    free(m); 
    *lp = l; 
}

static void
print_time(double const seconds)
{
    printf("Search Time: %0.06fs\n", seconds);
}

static void
print_numbers(
        char const * const filename,
        int const n,
        float const * const numbers)
{
    int i;
    FILE * fout;

    /* open file */
    if(NULL == (fout = fopen(filename, "w"))) {
        fprintf(stderr, "error opening '%s'\n", filename);
        abort();
    }

    /* write numbers to fout */
    for(i=0; i<n; ++i) {
        fprintf(fout, "%10.4f\n", numbers[i]);
    }

    fclose(fout);
}

int
main(int argc, char ** argv)
{
    int n, npe, rank;
    double ts, te;
    float * a, * l;
    int * offset = NULL, *nlocal = NULL; 
    if(argc < 4){
        printf("Invalid number of arguments.\nUsage: dijkstra <graph> <source> <output_file>.\n");
        return EXIT_FAILURE;
    }

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &npe);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == MAIN_PROCESS) {
        read_file_and_send(argv[1], &n, &a, npe, &offset, &nlocal);
    } else {
        recv_values_from_master(&n, &a, npe, &offset, &nlocal);
    }
    
    l = malloc(n*sizeof(*l));
    assert(l);
    ts = MPI_Wtime();
    dijkstra(atoi(argv[2]), n, a, &l, rank, offset, nlocal, npe);
    te = MPI_Wtime();
    if (rank == MAIN_PROCESS) {
        print_time(te - ts);
        print_numbers(argv[3], n, l);
    }
    
    free(a);
    free(l);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
