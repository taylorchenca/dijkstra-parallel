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

#define MASTER 0
#define N_TAG 0
#define A_TAG 1

static int get_chunk_size(int npe, int n, int rank) {
    if (npe > n) {
        if (rank > (npe - 1)) {
            return 0; 
        } else {
            return 1; 
        }
    }
    int chunk_size = n / npe; 
    int remainder = n % npe; 
    if (rank == npe - 1) return chunk_size + remainder;
    else return chunk_size; 
}

static int get_offset(int npe, int n, int rank) {
    int chunk_size = n / npe; 
    return rank * chunk_size; 
}

static int get_source_node(int npe, int n, int s) {
    int chunk_size = n / npe; 
    int i; 
    for (i = 0; i < npe; i++) {
        if (s >= i * chunk_size && s < (i + 1) * chunk_size) {
            return i; 
        }
    }
    return i; 
}

static void
read_file_and_send(
        char const * const filename,
        int * const np,
        float ** const ap,
        int npe
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

    int chunk_size;
    
    chunk_size = get_chunk_size(npe, n, 0);
    a = malloc(n * chunk_size * sizeof(float));
    for (j = 0; j < chunk_size * n; j++) {
        ret = fscanf(fp, "%f", &a[j]);
        assert(1 == ret);
    }
    *ap = a;

    for (i = 1; i < npe; i++) {
        chunk_size = get_chunk_size(npe, n, i);
        a = malloc(n * chunk_size * sizeof(float));
        MPI_Send(&n, 1, MPI_INTEGER, i, N_TAG, MPI_COMM_WORLD);

        for (j = 0; j < chunk_size * n; j++) {
            ret = fscanf(fp, "%f", &a[j]);
            assert(1 == ret);
        }
        MPI_Send(a, (chunk_size * n), MPI_FLOAT, i, A_TAG, MPI_COMM_WORLD);
        
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
        int rank
)
{
    int n;
    float *a = NULL; 
    MPI_Recv(&n, 1, MPI_INTEGER, 0, N_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    int chunk_size = get_chunk_size(npe, n, rank); 
    a = malloc(n * chunk_size * sizeof(float)); 
    MPI_Recv(a, (chunk_size * n), MPI_FLOAT, MASTER, A_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
    *ap = a; 
    *np = n; 
    MPI_Barrier(MPI_COMM_WORLD); 
}

static void
dijkstra(
        int const s, // Source
        int const n, // Number of vertices
        float const * const a, // Adjacency matrix
        float **const lp, // Result
        int rank, 
        int npe,  /* Number of compute nodes */ 
        int source_node
)
{
    int i, j;
    struct float_int {
        float l;
        int u;
    } min;
    
    char * m; /* m[i] == 1: i is visited */
    float *global_l = NULL;
    float *local_l = NULL;

    m = calloc(n, sizeof(*m));
    assert(m);
    global_l = malloc(n * sizeof(float)); 
    assert(global_l); 
    local_l = malloc(n * sizeof(float)); 
    assert(local_l);  

    int chunk_size = get_chunk_size(npe, n, rank); 
    int offset = get_offset(npe, n, rank); 

    if (rank == source_node) {
        int source_node_offset = get_offset(npe, n, source_node); 
        for (i = 0; i < n; i++) {
            /* Initialize distances to the distance from source to all other vertices */ 
            global_l[i] = a(s - source_node_offset, i); 
        }
    }
    
    MPI_Bcast(global_l, n, MPI_FLOAT, source_node, MPI_COMM_WORLD); /* Broadcast from source node to all others */ 
    MPI_Barrier(MPI_COMM_WORLD); 

    m[s] = 1; /* source vertex visited */ 
    min.u = -1; /* avoid compiler warning */

    for (i = 1; i < n; i++) {
        min.l = INFINITY;
        for (j = 0; j < n; j++) {
            if (!m[j] && global_l[j] < min.l) {
                min.l = global_l[j];
                min.u = j;
            }
            local_l[j] = global_l[j]; 
        }

        m[min.u] = 1; 
        for (j = 0; j < chunk_size; j++) {
            /* If a shorter path is found */ 
            if (!m[j + offset] && a(j, min.u) + min.l < local_l[j + offset]) {
                local_l[j + offset] = a(j, min.u) + min.l;
            }
        }
        /* Send out local_l for MIN reduction. Final result is received in global_l */ 
        MPI_Allreduce(local_l, global_l, n, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD); 
        MPI_Barrier(MPI_COMM_WORLD); 
    }
    free(m); 
    *lp = global_l; 
}

static void
print_time(double const seconds)
{
    printf("%0.06fs\n", seconds);
}

static void
print_log(        
    char const * const filename,
    int const npe,
    int const n,
    int const source,
    double const seconds
)
{
    FILE * fout;

    /* open file */
    if(NULL == (fout = fopen(filename, "a"))) {
        fprintf(stderr, "error opening '%s'\n", filename);
        abort();
    }

    fprintf(fout, "%d, %d, %d, %0.06fs\n", npe, n, source, seconds);
    fclose(fout);
    
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
    if(argc < 4){
        printf("Invalid number of arguments.\nUsage: dijkstra <graph> <source> <output_file> <log file>.\n");
        return EXIT_FAILURE;
    }

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &npe);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == MASTER) {
        read_file_and_send(argv[1], &n, &a, npe);
    } else {
        recv_values_from_master(&n, &a, npe, rank);
    }

    int s = atoi(argv[2]); 
    int source_node = get_source_node(npe, n, s); 
    l = malloc(n*sizeof(*l));
    assert(l);

    ts = MPI_Wtime();
    dijkstra(s, n, a, &l, rank, npe, source_node);
    te = MPI_Wtime();
    if (rank == MASTER) {
        print_log(argv[4], npe, n, s, te - ts); 
        print_time(te - ts);
        print_numbers(argv[3], n, l);
    }
    
    free(a);
    free(l);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
