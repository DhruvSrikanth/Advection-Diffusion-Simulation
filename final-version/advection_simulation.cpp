#include <stdio.h>
#include <iostream>
using namespace std;
#include <fstream>

#include <assert.h>

#include <string>

#include <math.h>

#include <mpi.h>
#include<omp.h>

// Dimensions of the distributed processor grid
// We have a 2D grid of processors
#define DIMENSION 2

void write_to_file(double** const& C_n, string const& filename, int const& NT, int const& N) { 
    // Allocate memory for the file
    ofstream file;
    file.open(filename);

    // Write the timestep to the file
    file << "[" << NT << "], ";

    // Write the data to the file
    file << "[";
    for (int i = 0; i < N; i++) {
        if (i == 0) {
            file << "[";
        } 
        else {
            file << ", [";
        }
        for (int j = 0; j < N; j++) {
            if (j == N - 1) {
                file << C_n[i][j];
            } 
            else {
                file << C_n[i][j] << ", ";
            }
        }
        file << "]";
    }
    file << "]";

    // Release the memory for the file
    file.close();
}

double** create_matrix(int const& N){
    // Allocate memory for the matrix
    double** mat = new double*[N];

    // Initialize the matrix
    for (int i = 0; i < N; i++) {
        mat[i] = new double[N];
        for (int j = 0; j < N; j++) {
            mat[i][j] = 0;
        }
    }
    return mat;
}

double** initial_gaussian(double** C_n, int const& N, double const& L) {
    int i;
    int j;
    // Using multiple threads, initialize the gaussian matrix
    #pragma omp parallel for default(none) private(i, j) shared(C_n, N, L) schedule(guided) 
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            C_n[i][j] = exp(-(pow((i-N/2)*L/N, 2) + pow((j-N/2)*L/N, 2))/(L*L/8));
        }
    }

    return C_n;
}

void apply_boundary_conditions(int const& N, double const& dt, double const& dx, double const& u, double const& v, double** C_n, double** C_n_1, double* const& left_ghost_cells, double* const& right_ghost_cells, double* const& up_ghost_cells, double* const& down_ghost_cells) {
    // Initialize the variables
    double C_n_up;
    double C_n_down;
    double C_n_left;
    double C_n_right;
    int i;
    int j;

    // Use multiple threads to apply the boundary conditions
    #pragma omp parallel for default(none) private(i, j, C_n_up, C_n_down, C_n_left, C_n_right) shared(C_n, C_n_1, N, dt, dx, u, v, left_ghost_cells, right_ghost_cells, up_ghost_cells, down_ghost_cells) schedule(guided)
    for (i = 0; i < N; i ++) {
        for (j = 0; j < N; j ++) {
            // Apply periodic boundary conditions
            // Middle rows
            if (i != 0 && i != N-1) {
                C_n_up = C_n[i-1][j];
                C_n_down = C_n[i+1][j];
                // Middle columns
                if (j != 0 && j != N-1) {
                    C_n_left = C_n[i][j-1];
                    C_n_right = C_n[i][j+1];
                }
                // Leftmost column
                else if (j == 0) {
                    //C_n_left = C_n[i][N-1];
                    C_n_left = left_ghost_cells[i];
                    C_n_right = C_n[i][j+1];
                } 
                // Rightmost column
                else {
                    C_n_left = C_n[i][j-1];
                    // C_n_right = C_n[i][0];
                    C_n_right = right_ghost_cells[i];
                }
            }
            // Top row
            else if (i == 0) {
                // C_n_up = C_n[N-1][j];
                C_n_up = up_ghost_cells[j];
                C_n_down = C_n[i+1][j];
                // Middle columns
                if (j != 0 && j != N-1) {
                    C_n_left = C_n[i][j-1];
                    C_n_right = C_n[i][j+1];
                }
                // Leftmost column
                else if (j == 0) {
                    // C_n_left = C_n[i][N-1];
                    C_n_left = left_ghost_cells[i];
                    C_n_right = C_n[i][j+1];
                } 
                // Rightmost column
                else {
                    C_n_left = C_n[i][j-1];
                    // C_n_right = C_n[i][0];
                    C_n_right = right_ghost_cells[i];
                }
            } 
            // Bottom row
            else {
                C_n_up = C_n[i-1][j];
                // C_n_down = C_n[0][j];
                C_n_down = down_ghost_cells[j];
                // Middle columns
                if (j != 0 && j != N-1) {
                    C_n_left = C_n[i][j-1];
                    C_n_right = C_n[i][j+1];
                }
                // Leftmost column
                else if (j == 0) {
                    // C_n_left = C_n[i][N-1];
                    C_n_left = left_ghost_cells[i];
                    C_n_right = C_n[i][j+1];
                } 
                // Rightmost column
                else {
                    C_n_left = C_n[i][j-1];
                    // C_n_right = C_n[i][0];
                    C_n_right = right_ghost_cells[i];
                }
            }
            
            // Compute next time step
            // Based on LAX formulation
            C_n_1[i][j] = 0.25*(C_n_up + C_n_down + C_n_left + C_n_right) - (dt/(2*dx))*(u*(C_n_down - C_n_up) + v*(C_n_right - C_n_left));
        }
    }
}

void advection_simulation(int const& N, int const& NT, double const& L, double const& T, double const& u, double const& v, int const& up, int const& down, int const& left, int const& right, MPI_Comm const& comm2d, int const& mype, int const& nprocs_per_dim, int const& nprocs) {
     
    // Initialize variables
    double** C_n = create_matrix(N);
    double** C_n_1 = create_matrix(N);
    double** global_output = create_matrix(N*nprocs_per_dim);

    double dx = L/N;
    double dt = T/NT;

    double best_grind_rate = 0.0;

    // Check Courant stability condition
    assert(dt<=dx/sqrt(2*(u*u + v*v)));

    // Initialize gaussian grid
    C_n = initial_gaussian(C_n, N, L);

    // Initialize columns to send
    double col_1[N];
    double col_n[N];

    for (int i = 0; i < N; i++) {
        col_1[i] = C_n[i][0];
        col_n[i] = C_n[i][N-1];
    }

    MPI_Status status;
    
    // Initialize ghost cells
    double up_ghost_cells[N];
    double down_ghost_cells[N];
    double left_ghost_cells[N];
    double right_ghost_cells[N];

    // Send and receive ghost cells
    MPI_Sendrecv(C_n[0], N, MPI_DOUBLE, up, 99, &down_ghost_cells, N, MPI_DOUBLE, down, MPI_ANY_TAG, comm2d, &status);
    MPI_Sendrecv(C_n[N-1], N, MPI_DOUBLE, down, 99, &up_ghost_cells, N, MPI_DOUBLE, up, MPI_ANY_TAG, comm2d, &status);
    MPI_Sendrecv(&col_1, N, MPI_DOUBLE, left, 99, &right_ghost_cells, N, MPI_DOUBLE, right, MPI_ANY_TAG, comm2d, &status);
    MPI_Sendrecv(&col_n, N, MPI_DOUBLE, right, 99, &left_ghost_cells, N, MPI_DOUBLE, left, MPI_ANY_TAG, comm2d, &status);
    
    // Apply BCs
    apply_boundary_conditions(N, dt, dx, u, v, C_n, C_n_1, left_ghost_cells, right_ghost_cells, up_ghost_cells, down_ghost_cells);
    
    // Swap references
    swap(C_n, C_n_1);

    // Write to file
    write_to_file(C_n, "./final-version/initial_gaussian.txt", nprocs_per_dim, N);

    for (int n = 0; n < NT; n++) {

        // Initialize columns to send
        for (int i = 0; i < N; i++) {
            col_1[i] = C_n[i][0];
            col_n[i] = C_n[i][N-1];
        }

        MPI_Status status;

        // Send and receive ghost cells
        MPI_Sendrecv(C_n[0], N, MPI_DOUBLE, up, 99, &down_ghost_cells, N, MPI_DOUBLE, down, MPI_ANY_TAG, comm2d, &status);
        MPI_Sendrecv(C_n[N-1], N, MPI_DOUBLE, down, 99, &up_ghost_cells, N, MPI_DOUBLE, up, MPI_ANY_TAG, comm2d, &status);
        MPI_Sendrecv(&col_1, N, MPI_DOUBLE, left, 99, &right_ghost_cells, N, MPI_DOUBLE, right, MPI_ANY_TAG, comm2d, &status);
        MPI_Sendrecv(&col_n, N, MPI_DOUBLE, right, 99, &left_ghost_cells, N, MPI_DOUBLE, left, MPI_ANY_TAG, comm2d, &status);

        // Start timer
	    double ts = MPI_Wtime();

        // Apply BCs
        apply_boundary_conditions(N, dt, dx, u, v, C_n, C_n_1, left_ghost_cells, right_ghost_cells, up_ghost_cells, down_ghost_cells);
        
        // Swap references
        swap(C_n, C_n_1);

        // Stop timer
	    double te = MPI_Wtime();
        
        // Output grind rate
        cout << "Iteration: " << n << " - Grind Rate: " << floor(1/(te-ts)) << " iter/sec" << endl;
        // Compute best grind rate
        best_grind_rate = max(best_grind_rate, floor(1/(te-ts)));

        // Write to file at the halfway point (time step = NT/2)
        if (n == (NT/2) - 1) {

            // Initialize local buffer for each processor
            double** local_buffer = create_matrix(N);

            // For processor 0
            if (mype == 0) {
                // Processor 0's contribution to the global output
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        global_output[i][j] = C_n[i][j];
                    }
                }

                // Every other processor's contribution to the global output
                for (int i = 1; i < nprocs; i++) {
                    MPI_Recv(*local_buffer, pow(N, 2), MPI_DOUBLE, i, 0, comm2d, MPI_STATUS_IGNORE);
                    for (int j = N * (i % nprocs_per_dim); j < N * ((i % nprocs_per_dim)+1); j++) {
                        for (int k = N * floor(i / nprocs_per_dim); k < N * (floor(i / nprocs_per_dim) + 1); k++) {
                            global_output[j][k] = local_buffer[j - (int(N) * (i % nprocs_per_dim))][k - (int(N * floor(i / nprocs_per_dim)))];
                        }
                    }
                }

                // Write to file
                write_to_file(global_output, "./final-version/simulation_NTby2_timesteps.txt", nprocs_per_dim, N);
            }
            // For other processors
            else {
                MPI_Send(*C_n, pow(N, 2), MPI_DOUBLE, 0, 0, comm2d);
            }
        } 

        // Write to file upon completion of simulation (time step = NT)
        else if (n == NT - 1) {

            // Initialize local buffer for each processor
            double** local_buffer = create_matrix(N);

            // For processor 0
            if (mype == 0) {
                // Processor 0's contribution to the global output
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        global_output[i][j] = C_n[i][j];
                    }
                }

                // Every other processor's contribution to the global output
                for (int i = 1; i < nprocs; i++) {
                    MPI_Recv(*local_buffer, pow(N, 2), MPI_DOUBLE, i, 0, comm2d, MPI_STATUS_IGNORE);
                    for (int j = N * (i % nprocs_per_dim); j < N * ((i % nprocs_per_dim)+1); j++) {
                        for (int k = N * floor(i / nprocs_per_dim); k < N * (floor(i / nprocs_per_dim) + 1); k++) {
                            global_output[j][k] = local_buffer[j - (int(N) * (i % nprocs_per_dim))][k - (int(N * floor(i / nprocs_per_dim)))];
                        }
                    }
                }

                // Write to file
                write_to_file(global_output, "./final-version/simulation_NT_timesteps.txt", nprocs_per_dim, N);
            }
            // For other processors 
            else {
                MPI_Send(*C_n, pow(N, 2), MPI_DOUBLE, 0, 0, comm2d);
            }
        }
    }

    // Best grind rate
    cout << "Best grind rate: " << best_grind_rate << " iter/sec" << endl;
}

int main(int argc, char** argv) {
    // MPI initialization
	int mype, nprocs;

    MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    // MPI Cartesian Grid Creation
	int dims[DIMENSION], periodic[DIMENSION], coords[DIMENSION];
	int nprocs_per_dim = floor(sqrt(nprocs));
	MPI_Comm comm2d;
	dims[0] = nprocs_per_dim;
    dims[1] = nprocs_per_dim;
	periodic[0] = 1;
    periodic[1] = 1;

    // Create Cartesian Communicator
	MPI_Cart_create(MPI_COMM_WORLD, DIMENSION, dims, periodic, 1, &comm2d);

	// Extract this MPI rank's N-dimensional coordinates from its place in the MPI Cartesian grid
	MPI_Cart_coords(comm2d, mype, DIMENSION, coords);

	// Determine 2D neighbor ranks for this MPI rank
	int left, right;
	MPI_Cart_shift(comm2d, 1, 1, &left, &right);
    
    int down, up;
	MPI_Cart_shift(comm2d, 0, 1, &down, &up);


    // Initialize variables
    int N_glb = stoi(argv[1]);
    int N_loc = N_glb/nprocs_per_dim;

    int NT = stoi(argv[2]);
    double L = stod(argv[3]);
    double T = stod(argv[4]);
    double u = stod(argv[5]);
    double v = stod(argv[6]);

    int num_threads;

    cout << "Simulation Parameters:" << endl;
    cout << "N = " << N_glb << endl;
    cout << "NT = " << NT << endl;
    cout << "L = " << L << endl;
    cout << "T = " << T << endl;
    cout << "u = " << u << endl;
    cout << "v = " << v << "\n" << endl;
    
    // Esimate memory usage
    cout << "Estimated memeory usage = " << N_loc*N_loc*sizeof(double)/1e6 << " MB" << "\n" << endl;

    // Set number of threads
    num_threads = 6;
    cout << "Number of threads = " << num_threads << "\n" << endl;
    omp_set_num_threads(num_threads);

    // Start timer
    double t1 = omp_get_wtime();

    // Run simulation
    cout << "Simulating..." << endl;
    advection_simulation(N_loc, NT, L, T, u, v, up, down, left, right, comm2d, mype, nprocs_per_dim, nprocs);
    cout << "Simulation Complete!" << endl;

    // Stop timer
    double t2 = omp_get_wtime();

    // MPI finalization
    MPI_Finalize();

    // Print time taken
    cout << "Time taken = " << t2-t1 << " seconds" << "\n" << endl;
}