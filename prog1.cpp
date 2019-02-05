#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
using namespace std;

double randfrom(double mini, double maxi) {
  double range = maxi - mini;
  double div = RAND_MAX / range;
  return mini + (rand() / div);
}

int dboard(int N) {
  MPI_Comm comm = MPI_COMM_WORLD;
  int p;
  int darts_inside_square = 0;
  MPI_Comm_size(comm, &p);
  int total_darts = floor(N / p);

  for (int i = 0; i < total_darts; i++) {
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

      /* code */
      if (argc != 3) {
        return 0;
      }

      int N, R;

      if (rank == 0) {
        N = atoi(argv[1]);
        R = atoi(argv[2]);
      }

      MPI_Bcast(&N, 1, MPI_INT, 0, comm);
      MPI_Bcast(&R, 1, MPI_INT, 0, comm);

      double sum_estimates = 0.0;
      double runtime = 0;
      int M;
      for (int i = 0; i < R; i++) {
        double t0 = MPI_Wtime();
        srand(rank + 1);
        int m = dboard(N);
        MPI_Reduce(&m, &M, 1, MPI_INT, MPI_SUM, 0, comm);
        double t1 = MPI_Wtime();
        runtime += (t1 - t0);

        if (rank == 0) {
          double pi_estimate = ((double) N / M) * 2.0;
          sum_estimates += pi_estimate;
        }

      }

      double final_estimate = sum_estimates / (double) R;
      if (rank == 0) {
        std::cout << "N = " << N << "  R = " << R << " P = " << p << " PI = " << final_estimate << "\n";
        std::cout << "Time = " << runtime << "\n";
      }
      
      // if (rank == 0) {
      //   printf("num arguments = %d, N = %d, R = %d, final estimate = %f, rank = %d, runtime = %f \n", argc, N, R, final_estimate, rank, runtime);
      // }
    // finalize MPI
    MPI_Finalize();
    return 0;
}
