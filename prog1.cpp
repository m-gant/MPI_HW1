#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
using namespace std;

//creates a random double inbetween a min and max value
double randfrom(double min, double max) {
  double range = max - min;
  double div = RAND_MAX / range;
  return min + (rand() / div);
}

/*returns the amount of darts that land inside of the largest square
that lies within a circle*/
int dboard(int N) {
  MPI_Comm comm = MPI_COMM_WORLD;
  int p;
  int darts_inside_square = 0;
  MPI_Comm_size(comm, &p);

  //calculating the amount of darts to simulate with one processor
  int total_darts = floor(N / p);

  for (int i = 0; i < total_darts; i++) {
      //using a circle with radius 2
      double alpha = randfrom(0.0, 4.0);
      double r = sqrt(alpha);
      double theta = randfrom(0.0, 360.0);
      double x = r * cos(theta * M_PI / 180.0);
      double y = r * sin(theta * M_PI / 180.0);
      if ((abs(x) <= (2.0 / (sqrt(2)))) && (abs(y) <= (2.0 / sqrt(2)))) {
        darts_inside_square++;
      }
  }
  return darts_inside_square;
}


int main(int argc,char*argv[]){
  // set up MPI
  MPI_Init(&argc,&argv);
  // get communicator size and my rank
      MPI_Comm comm = MPI_COMM_WORLD;
      int p, rank;
      MPI_Comm_size(comm,&p);
      MPI_Comm_rank(comm,&rank);

      //checking if there are the right amount of input variables
      if (argc != 3) {
        return 0;
      }

      int N, R;

      if (rank == 0) {
        N = atoi(argv[1]);
        R = atoi(argv[2]);
      }

      //broadcasting the input variables to all of the non-root processors
      MPI_Bcast(&N, 1, MPI_INT, 0, comm);
      MPI_Bcast(&R, 1, MPI_INT, 0, comm);

      double sum_estimates = 0.0;
      double runtime = 0;
      int M;

      for (int i = 0; i < R; i++) {
        double t0 = MPI_Wtime();
        //seeding random according to the rank of the processor.
        srand(rank + 1);
        int m = dboard(N);
        MPI_Reduce(&m, &M, 1, MPI_INT, MPI_SUM, 0, comm);
        double t1 = MPI_Wtime();
        runtime += (t1 - t0);

        if (rank == 0) {
          //adding up all of the darts and calculating pi for this iteration.
          double pi_estimate = ((double) N / M) * 2.0;
          sum_estimates += pi_estimate;
        }

      }

      //calculating the average of all R iterations of our dartboard rounds.
      double final_estimate = sum_estimates / (double) R;

      //outpuing the results of the program
      if (rank == 0) {
        std::cout << "N = " << N << "  R = " << R << " P = " << p << " PI = " << final_estimate << "\n";
        std::cout << "Time = " << runtime << "\n";
      }


    // finalize MPI
    MPI_Finalize();
    return 0;
}
