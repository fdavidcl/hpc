#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

struct timespec diff(struct timespec start, struct timespec end) {
  struct timespec temp;
  if ((end.tv_nsec - start.tv_nsec) < 0) {
    temp.tv_sec = end.tv_sec - start.tv_sec - 1;
    temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec - start.tv_sec;
    temp.tv_nsec = end.tv_nsec - start.tv_nsec;
  }
  return temp;
}

double approximate_pi(double width, long ini, long end, double offset) {
  double sum = 0;

  long i;
  //#pragma omp for
  for (i = ini; i < end; ++i) {
    double x = (i + offset) * width;
    sum += 1 / (1.0 + x * x);
  }

  return 4.0 * sum * width;
}
  
int main(int argc, char **argv) {
  const double PI = 3.141592653589793238462643383279502884197169399375105820974; // wolframalpha rules

  // Approximation methods
  const double
    EXCESS = 1,
    DEFECT = 0,
    MIDDLE = 0.5;
  
  if (argc < 2) {
    printf("Mu mal\n");
    return 1;
  }
  
  long intervals = atol(argv[1]);
  double width = 1.0 / intervals;
  double
    errore = PI,
    errord = PI,
    errorm = PI;

  /*** Initialize MPI ***/
  const int MASTER_RANK = 0;
  int size, rank;
  
  MPI_Status status;

  struct timespec timprocini, timprocend;
  clock_gettime(CLOCK_REALTIME, &timprocini);
  
  MPI_Init(&argc, &argv);

  clock_gettime(CLOCK_REALTIME, &timprocend);
  struct timespec timproc = diff(timprocini, timprocend),
    timcalc, timcomms;
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /*
  if (size < 2) {
    printf("Necesito al menos 2 procesos");
    MPI_Finalize();
    return 1;
    }*/
  
  long
    /*** Calculate how many intervals per process ***/
    ints_per_proc = intervals / size,
    /* The first interval for this process */
    ini = rank * ints_per_proc,
    /* The next to the last interval for this process */
    end = (rank < size - 1) ? (rank + 1) * ints_per_proc : intervals;

  struct timespec tcalcini, tcalcend;
  
  double partials[3];

  if (rank == MASTER_RANK)
    clock_gettime(CLOCK_REALTIME, &tcalcini);

  partials[0] = approximate_pi(width, ini, end, EXCESS);
  partials[1] = approximate_pi(width, ini, end, DEFECT);
  partials[2] = approximate_pi(width, ini, end, MIDDLE);
  
  if (rank == MASTER_RANK) {
    clock_gettime(CLOCK_REALTIME, &tcalcend);
    timcalc = diff(tcalcini, tcalcend);
  }

  //printf("I'm rank %d and I calculated %f\n", rank, partials[2]);

  if (rank == MASTER_RANK) {
    errore -= partials[0];
    errord -= partials[1];
    errorm -= partials[2];

    struct timespec tcommini, tcommend;
    clock_gettime(CLOCK_REALTIME, &tcommini);
    
    /* Receive all partial sums */
    int comp;
    for (comp = 1; comp < size; comp++) {
      MPI_Recv(partials, 3, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      errore -= partials[0];
      errord -= partials[1];
      errorm -= partials[2];
    }

    clock_gettime(CLOCK_REALTIME, &tcommend);
    timcomms = diff(tcommini, tcommend);
  } else {
    MPI_Send(partials, 3, MPI_DOUBLE, MASTER_RANK, 0, MPI_COMM_WORLD);
  }

  clock_gettime(CLOCK_REALTIME, &timprocini);
  
  MPI_Finalize();
  
  clock_gettime(CLOCK_REALTIME, &timprocend);
  struct timespec timproc2 = diff(timprocini, timprocend);

  if (rank == MASTER_RANK) {
    /***
        Formato:
        número de procesos, número de intervalos, error por exceso, error por defecto, error en punto medio, tiempo de cálculo, tiempo de comunicaciones, tiempo de creación de procesos 
    ***/
    printf("%d,%li,%.20f,%.20f,%.20f,%i.%09li,%i.%09li,%i.%09li,%i.%09li",
           size,
           intervals,
           errore,
           errord,
           errorm,
           timcalc.tv_sec,
           timcalc.tv_nsec,
           timcomms.tv_sec,
           timcomms.tv_nsec,
           timproc.tv_sec,
           timproc.tv_nsec,
           timproc2.tv_sec,
           timproc2.tv_nsec
           );
  }

  return 0;
}

