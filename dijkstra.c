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

#define ROWMJR(R,C,NR,NC) (R*NC+C)
#define COLMJR(R,C,NR,NC) (C*NR+R)
/* define access directions for matrices */
#define a(R,C) a[ROWMJR(R,C,ln,n)]
#define b(R,C) b[ROWMJR(R,C,nn,n)]

static void
load(
  char const * const filename,
  int * const np,
  float ** const ap
)
{
  int i, j, n, ret;
  FILE * fp=NULL;
  float * a;

  /* open the file */
  fp = fopen(filename, "r");
  assert(fp);

  /* get the number of nodes in the graph */
  ret = fscanf(fp, "%d", &n);
  assert(1 == ret);

  /* allocate memory for local values */
  a = malloc(n*n*sizeof(*a));
  assert(a);

  /* read in roots local values */
  for (i=0; i<n; ++i) {
    for (j=0; j<n; ++j) {
      ret = fscanf(fp, "%f", &a(i,j));
      assert(1 == ret);
    }
  }

  /* close file */
  ret = fclose(fp);
  assert(!ret);

  /* record output values */
  *np = n;
  *ap = a;
}

static void
dijkstra(
  int const s, // Source
  int const n, // Number of nodes
  float const * const a, // Adjacency matrix
  float ** const lp // Result
)
{
  int i, j;
  struct float_int {
    float l;
    int u;
  } min;
  /* m keeps track of visited. m[v] == 1 => v has been visited and no longer in the queue */
  char * m;
  /* l[v]: distance of source to v */
  float * l;

  m = calloc(n, sizeof(*m));
  assert(m);

  l = malloc(n*sizeof(*l));
  assert(l);

  for (i=0; i<n; ++i) {
    l[i] = a(i,s); /* Initialize l[v] */
  }


  m[s] = 1; /* Source is visited. */
  min.u = -1; /* avoid compiler warning */

  for (i=1; i<n; ++i) {
    min.l = INFINITY;

    /* find local minimum
     * Find j such that j is the smallest among l and not yet visited */
    for (j=0; j<n; ++j) {
      if (!m[j] && l[j] < min.l) {
        min.l = l[j];
        min.u = j;
      }
    }

    m[min.u] = 1; /* Mark the j found in the prior for loop as visited */

    /*
     * Find alternative routes to minimize l[j]
     */
    for (j=0; j<n; ++j) {
      if (!m[j] && min.l+a(j,min.u) < l[j])
        l[j] = min.l+a(j,min.u);
    }
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
  int n;
  clock_t ts, te;
  float * a, * l;

  if(argc < 4){
     printf("Invalid number of arguments.\nUsage: dijkstra <graph> <source> <output_file>.\n");
     return EXIT_FAILURE;
  }


  load(argv[1], &n, &a);

  ts = clock();
  dijkstra(atoi(argv[2]), n, a, &l);
  te = clock();

  print_time((double)(te-ts)/CLOCKS_PER_SEC);
  print_numbers(argv[3], n, l);

  free(a);
  free(l);

  return EXIT_SUCCESS;
}
