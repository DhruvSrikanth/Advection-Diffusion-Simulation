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

void write_to_file(double** C_n, string const& filename, int const& NT, int const& N) { 
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

int contiguous_memory_index(int const& i, int const& j, int const& N) {
    // Calculate the index of the contiguous memory matrix
    return j + (N * i);
}

double* contiguous_memory_alloc(int const& N) {
    // Allocate contiguous memory for the matrix
    double* mat = new double[N * N];
    // Initialize the matrix
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            mat[contiguous_memory_index(i, j, N)] = 0.0;
        }
    }
    return mat;
}

void initial_gaussian(double** C_n, int const& N, double const& L, int const& N_glob, int const& mype, int const& nprocs_per_dim) {

    // Initialize variables
    int i;
    int j;
    int glob_i;
    int glob_j; 

    // Calculate the global start (i,j) with respect to the current processor
    int glob_i_start = (nprocs_per_dim - 1 - floor(mype / nprocs_per_dim)) * N;
    int glob_j_start = (mype % nprocs_per_dim) * N;

    // Using multiple threads, initialize the gaussian matrix
    #pragma omp parallel for default(none) private(i, j, glob_i, glob_j) shared(C_n, N, L, mype, nprocs_per_dim, N_glob, glob_i_start, glob_j_start) schedule(guided) 
    for (i = 0; i < N; i++) {
        // Calculate the global i
        glob_i =  glob_i_start + i;
        for (j = 0; j < N; j++) {
            // Calculate the global j
            glob_j =  glob_j_start + j;
            // Calculate the gaussian value with respect to the global coordinates
            C_n[i][j] = exp(-(pow((glob_j-(N_glob/2))*L/N_glob, 2) + pow((glob_i-(N_glob/2))*L/N_glob, 2))/(L*L/8));
        }
    }

    // Helpful for debugging
    // Save each processor's data to a file
    // char filename[100];
    // cout << "Writing output to file..." << endl;
    // int dummy_var = sprintf(filename, "./final-version/mype_%d.txt", mype);
    // write_to_file(C_n, filename, 0, N);
}

void apply_boundary_conditions(int const& N, double const& dt, double const& dx, double const& u, double const& v, double** C_n, double** C_n_1, double* const& left_ghost_cells, double* const& right_ghost_cells, double* const& up_ghost_cells, double* const& down_ghost_cells, string const& scheme) {
    // Initialize the variables accounting for 1D and 2D ghost cells
    double C_n_up;
    double C_n_up2;
    double C_n_down;
    double C_n_down2;
    double C_n_left;
    double C_n_left2;
    double C_n_right;
    double C_n_right2;
    int i;
    int j;


    // Use multiple threads to apply the boundary conditions
    #pragma omp parallel for default(none) private(i, j, C_n_up, C_n_down, C_n_left, C_n_right, C_n_up2, C_n_down2, C_n_left2, C_n_right2) shared(C_n, C_n_1, N, dt, dx, u, v, left_ghost_cells, right_ghost_cells, up_ghost_cells, down_ghost_cells, scheme) schedule(guided)
    for (i = 0; i < N; i ++) {
        for (j = 0; j < N; j ++) {
            // Apply periodic boundary conditions
            if (scheme == "Second-Order-Upwind") {
                // Top row
                if (i == 0) {
                    C_n_up = up_ghost_cells[j];
                    C_n_up2 = up_ghost_cells[j + N];
                    C_n_down = C_n[i+1][j];
                    C_n_down2 = C_n[i+2][j];
                    // Left column
                    if (j == 0) {
                        C_n_left = left_ghost_cells[i];
                        C_n_left2 = left_ghost_cells[i + N];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                    // Right column
                    else if (j == N - 1) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = right_ghost_cells[i];
                        C_n_right2 = right_ghost_cells[i + N];
                    }
                    // Left + 1 column
                    else if (j == 1) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = left_ghost_cells[i];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                    // Right - 1 column
                    else if (j == N - 2) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = right_ghost_cells[i];
                    }
                    // Middle columns
                    else {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                }
                // Bottom row
                else if (i == N - 1) {
                    C_n_up = C_n[i-1][j];
                    C_n_up2 = C_n[i-2][j];
                    C_n_down = down_ghost_cells[j];
                    C_n_down2 = down_ghost_cells[j + N];
                    // Left column
                    if (j == 0) {
                        C_n_left = left_ghost_cells[i];
                        C_n_left2 = left_ghost_cells[i + N];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                    // Right column
                    else if (j == N - 1) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = right_ghost_cells[i];
                        C_n_right2 = right_ghost_cells[i + N];
                    }
                    // Left + 1 column
                    else if (j == 1) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = left_ghost_cells[i];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                    // Right - 1 column
                    else if (j == N - 2) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = right_ghost_cells[i];
                    }
                    // Middle columns
                    else {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }

                }
                // Top - 1 row
                else if (i == 1) {
                    C_n_up = C_n[i-1][j];
                    C_n_up2 = up_ghost_cells[j];
                    C_n_down = C_n[i+1][j];
                    C_n_down2 = C_n[i+2][j];
                    // Left column
                    if (j == 0) {
                        C_n_left = left_ghost_cells[i];
                        C_n_left2 = left_ghost_cells[i + N];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                    // Right column
                    else if (j == N - 1) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = right_ghost_cells[i];
                        C_n_right2 = right_ghost_cells[i + N];
                    }
                    // Left + 1 column
                    else if (j == 1) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = left_ghost_cells[i];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                    // Right - 1 column
                    else if (j == N - 2) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = right_ghost_cells[i];
                    }
                    // Middle columns
                    else {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                }
                // Bottom + 1 row
                else if (i == N - 2) {
                    C_n_up = C_n[i-1][j];
                    C_n_up2 = C_n[i-2][j];
                    C_n_down = C_n[i+1][j];
                    C_n_down2 = down_ghost_cells[j];
                    // Left column
                    if (j == 0) {
                        C_n_left = left_ghost_cells[i];
                        C_n_left2 = left_ghost_cells[i + N];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                    // Right column
                    else if (j == N - 1) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = right_ghost_cells[i];
                        C_n_right2 = right_ghost_cells[i + N];
                    }
                    // Left + 1 column
                    else if (j == 1) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = left_ghost_cells[i];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                    // Right - 1 column
                    else if (j == N - 2) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = right_ghost_cells[i];
                    }
                    // Middle columns
                    else {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                }
                // Middle rows
                else {
                    C_n_up = C_n[i-1][j];
                    C_n_up2 = C_n[i-2][j];
                    C_n_down = C_n[i+1][j];
                    C_n_down2 = C_n[i+2][j];
                    // Left column
                    if (j == 0) {
                        C_n_left = left_ghost_cells[i];
                        C_n_left2 = left_ghost_cells[i + N];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                    // Right column
                    else if (j == N - 1) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = right_ghost_cells[i];
                        C_n_right2 = right_ghost_cells[i + N];
                    }
                    // Left + 1 column
                    else if (j == 1) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = left_ghost_cells[i];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                    // Right - 1 column
                    else if (j == N - 2) {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = right_ghost_cells[i];
                    }
                    // Middle columns
                    else {
                        C_n_left = C_n[i][j-1];
                        C_n_left2 = C_n[i][j-2];
                        C_n_right = C_n[i][j+1];
                        C_n_right2 = C_n[i][j+2];
                    }
                }
            }
            else {
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
            }
            
            // Average (Smoothen) gaussian value
            C_n[i][j] = 0.25*(C_n_up + C_n_down + C_n_left + C_n_right);

            // Based on LAX formulation
            if (scheme == "LAX") {
                C_n_1[i][j] = C_n[i][j] - (dt/(2*dx))*(u*(C_n_down - C_n_up) + v*(C_n_right - C_n_left));
            }
            // Based on First Order Upwind Scheme
            else if (scheme == "First-Order-Upwind") {
                if (u*u + v*v > 0) { 
                    C_n_1[i][j] = C_n[i][j] - (dt/dx)*(u*(C_n[i][j] - C_n_up) + v*(C_n[i][j] - C_n_left));
                }
                else if (u*u + v*v < 0) {
                    C_n_1[i][j] = C_n[i][j] - (dt/dx)*(u*(C_n_down - C_n[i][j]) + v*(C_n_right - C_n[i][j]));
                }
            }
            // Based on Second Order Upwind Scheme
            else if (scheme == "Second-Order-Upwind") {
                if (u*u + v*v > 0) {
                    C_n_1[i][j] = C_n[i][j] - (dt/(2*dx))*(u*(3*C_n[i][j] - 4*C_n_up + C_n_up2) + v*(3*C_n[i][j] - 4*C_n_left + C_n_left2));
                }
                else if (u*u + v*v < 0) {
                    C_n_1[i][j] = C_n[i][j] - (dt/(2*dx))*(u*(4*C_n_down - C_n_down2 - 3*C_n[i][j]) + v*(4*C_n_right - 3*C_n[i][j] - C_n_right2));
                }
            }
        }
    }
}

void advection_simulation(int const& N, int const& NT, double const& L, double const& T, double const& u, double const& v, int const& up, int const& down, int const& left, int const& right, MPI_Comm const& comm2d, int const& mype, int const& nprocs_per_dim, int const& nprocs, string const& scheme) {
     
    // Initialize variables
    double** C_n = create_matrix(N);
    double** C_n_1 = create_matrix(N);
    int N_glob = N*nprocs_per_dim;

    // Only for processor 0 (master) to combine output from all processors
    double** global_output = create_matrix(N_glob);

    double dx = L/N_glob;
    double dt = T/NT;

    double best_grind_rate = 0.0;

    // Check Courant stability condition
    assert(dt<=dx/sqrt(2*(u*u + v*v)));

    // Initialize gaussian grid
    initial_gaussian(C_n, N, L, N_glob, mype, nprocs_per_dim);

    // Initialize data to send
    if (scheme == "Second-Order-Upwind") {
        double col_1[2*N];
        double col_n[2*N];
        double row_1[2*N];
        double row_n[2*N];
        int i;

        // Use multiple threads to compute data to send
        #pragma omp parallel for default(none) private(i) shared(row_1, row_n, col_1, col_n, C_n, N) schedule(guided)  
        for (i = 0; i < N; i++) {
            // left 2N columns
            col_1[i] = C_n[i][0];
            col_1[i+N] = C_n[i][1];
            
            // right 2N columns
            col_n[i] = C_n[i][N-1];
            col_n[i+N] = C_n[i][N-2];

            // top 2N rows
            row_1[i] = C_n[0][i];
            row_1[i+N] = C_n[1][i];

            // bottom 2N rows
            row_n[i] = C_n[N-1][i];
            row_n[i+N] = C_n[N-2][i];
        }

        double up_ghost_cells[2*N];
        double down_ghost_cells[2*N];
        double left_ghost_cells[2*N];
        double right_ghost_cells[2*N];

        MPI_Status status;

        // Send and receive ghost cells
        MPI_Sendrecv(&row_1, 2*N, MPI_DOUBLE, up, 99, &down_ghost_cells, 2*N, MPI_DOUBLE, down, MPI_ANY_TAG, comm2d, &status);
        MPI_Sendrecv(&row_n, 2*N, MPI_DOUBLE, down, 99, &up_ghost_cells, 2*N, MPI_DOUBLE, up, MPI_ANY_TAG, comm2d, &status);
        MPI_Sendrecv(&col_1, 2*N, MPI_DOUBLE, left, 99, &right_ghost_cells, 2*N, MPI_DOUBLE, right, MPI_ANY_TAG, comm2d, &status);
        MPI_Sendrecv(&col_n, 2*N, MPI_DOUBLE, right, 99, &left_ghost_cells, 2*N, MPI_DOUBLE, left, MPI_ANY_TAG, comm2d, &status);

        // Apply BCs
        apply_boundary_conditions(N, dt, dx, u, v, C_n, C_n_1, left_ghost_cells, right_ghost_cells, up_ghost_cells, down_ghost_cells, scheme);

    }
    else {
        double col_1[N];
        double col_n[N];
        int i;

        // Use multiple threads to compute data to send
        #pragma omp parallel for default(none) private(i) shared(col_1, col_n, C_n, N) schedule(guided)  
        for (i = 0; i < N; i++) {
            col_1[i] = C_n[i][0];
            col_n[i] = C_n[i][N-1];
        }
        
        double up_ghost_cells[N];
        double down_ghost_cells[N];
        double left_ghost_cells[N];
        double right_ghost_cells[N];

        MPI_Status status;

        // Send and receive ghost cells
        MPI_Sendrecv(C_n[0], N, MPI_DOUBLE, up, 99, &down_ghost_cells, N, MPI_DOUBLE, down, MPI_ANY_TAG, comm2d, &status);
        MPI_Sendrecv(C_n[N-1], N, MPI_DOUBLE, down, 99, &up_ghost_cells, N, MPI_DOUBLE, up, MPI_ANY_TAG, comm2d, &status);
        MPI_Sendrecv(&col_1, N, MPI_DOUBLE, left, 99, &right_ghost_cells, N, MPI_DOUBLE, right, MPI_ANY_TAG, comm2d, &status);
        MPI_Sendrecv(&col_n, N, MPI_DOUBLE, right, 99, &left_ghost_cells, N, MPI_DOUBLE, left, MPI_ANY_TAG, comm2d, &status);

        // Apply BCs
        apply_boundary_conditions(N, dt, dx, u, v, C_n, C_n_1, left_ghost_cells, right_ghost_cells, up_ghost_cells, down_ghost_cells, scheme);

    }
    
    // Swap references
    swap(C_n, C_n_1);

    // Initialize local buffer for each processor
    double* contiguous_local_buffer = contiguous_memory_alloc(N);

    // For processor 0
    if (mype == 0) {
        // Processor 0's contribution to the global output
        int glob_i;
        int glob_j; 
        for (int i = 0; i < N; i++) {
            glob_i =  ((nprocs_per_dim - 1) * N) + i;
            for (int j = 0; j < N; j++) {
                glob_j = j;
                global_output[glob_i][glob_j] = C_n[i][j];
            }
        }

        // Every other processor's contribution to the global output
        for (int x = 1; x < nprocs; x++) {
            MPI_Recv(contiguous_local_buffer, N*N, MPI_DOUBLE, x, 0, comm2d, MPI_STATUS_IGNORE);
            for (int i = 0; i < N; i++) {
                glob_i =  ((nprocs_per_dim - 1 - floor(x / nprocs_per_dim)) * N) + i;
                for (int j = 0; j < N; j++) {
                    glob_j = ((x % nprocs_per_dim) * N) + j;
                    global_output[glob_i][glob_j] = contiguous_local_buffer[contiguous_memory_index(i, j, N)];
                }
            }
        }
        // Write to file
        write_to_file(global_output, "./final-version/initial_gaussian.txt", nprocs_per_dim, N_glob);
        
    }
    // For other processors
    else {
        // Initialize local output as a contiguous array to send
        double* contiguous_C_n = contiguous_memory_alloc(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                contiguous_C_n[contiguous_memory_index(i, j, N)] = C_n[i][j];
            }
        }
        MPI_Send(contiguous_C_n, N*N, MPI_DOUBLE, 0, 0, comm2d);
    }

    // BREAKPOINT for debugging purposes
    // return ;    

    // Run simulation for a certain number of timesteps
    for (int n = 0; n < NT; n++) {

        // Start timer
	    double ts = MPI_Wtime();
        
        // Initialize data to send
        if (scheme == "Second-Order-Upwind") {
            double col_1[2*N];
            double col_n[2*N];
            double row_1[2*N];
            double row_n[2*N];
            int i;

            // Use multiple threads to compute data to send
            #pragma omp parallel for default(none) private(i) shared(row_1, row_n, col_1, col_n, C_n, N) schedule(guided)  
            for (i = 0; i < N; i++) {
                // left 2N columns
                col_1[i] = C_n[i][0];
                col_1[i+N] = C_n[i][1];
                
                // right 2N columns
                col_n[i] = C_n[i][N-1];
                col_n[i+N] = C_n[i][N-2];

                // top 2N rows
                row_1[i] = C_n[0][i];
                row_1[i+N] = C_n[1][i];

                // bottom 2N rows
                row_n[i] = C_n[N-1][i];
                row_n[i+N] = C_n[N-2][i];
            }

            double up_ghost_cells[2*N];
            double down_ghost_cells[2*N];
            double left_ghost_cells[2*N];
            double right_ghost_cells[2*N];

            MPI_Status status;

            // Send and receive ghost cells
            MPI_Sendrecv(&row_1, 2*N, MPI_DOUBLE, up, 99, &down_ghost_cells, 2*N, MPI_DOUBLE, down, MPI_ANY_TAG, comm2d, &status);
            MPI_Sendrecv(&row_n, 2*N, MPI_DOUBLE, down, 99, &up_ghost_cells, 2*N, MPI_DOUBLE, up, MPI_ANY_TAG, comm2d, &status);
            MPI_Sendrecv(&col_1, 2*N, MPI_DOUBLE, left, 99, &right_ghost_cells, 2*N, MPI_DOUBLE, right, MPI_ANY_TAG, comm2d, &status);
            MPI_Sendrecv(&col_n, 2*N, MPI_DOUBLE, right, 99, &left_ghost_cells, 2*N, MPI_DOUBLE, left, MPI_ANY_TAG, comm2d, &status);

            // Apply BCs
            apply_boundary_conditions(N, dt, dx, u, v, C_n, C_n_1, left_ghost_cells, right_ghost_cells, up_ghost_cells, down_ghost_cells, scheme);

        }
        else {
            double col_1[N];
            double col_n[N];
            int i;

            // Use multiple threads to compute data to send
            #pragma omp parallel for default(none) private(i) shared(col_1, col_n, C_n, N) schedule(guided)  
            for (i = 0; i < N; i++) {
                col_1[i] = C_n[i][0];
                col_n[i] = C_n[i][N-1];
            }
            
            double up_ghost_cells[N];
            double down_ghost_cells[N];
            double left_ghost_cells[N];
            double right_ghost_cells[N];

            MPI_Status status;

            // Send and receive ghost cells
            MPI_Sendrecv(C_n[0], N, MPI_DOUBLE, up, 99, &down_ghost_cells, N, MPI_DOUBLE, down, MPI_ANY_TAG, comm2d, &status);
            MPI_Sendrecv(C_n[N-1], N, MPI_DOUBLE, down, 99, &up_ghost_cells, N, MPI_DOUBLE, up, MPI_ANY_TAG, comm2d, &status);
            MPI_Sendrecv(&col_1, N, MPI_DOUBLE, left, 99, &right_ghost_cells, N, MPI_DOUBLE, right, MPI_ANY_TAG, comm2d, &status);
            MPI_Sendrecv(&col_n, N, MPI_DOUBLE, right, 99, &left_ghost_cells, N, MPI_DOUBLE, left, MPI_ANY_TAG, comm2d, &status);

            // Apply BCs
            apply_boundary_conditions(N, dt, dx, u, v, C_n, C_n_1, left_ghost_cells, right_ghost_cells, up_ghost_cells, down_ghost_cells, scheme);

        }
        
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
            double* contiguous_local_buffer = contiguous_memory_alloc(N);

            // For processor 0
            if (mype == 0) {
                // Processor 0's contribution to the global output
                int glob_i;
                int glob_j; 
                for (int i = 0; i < N; i++) {
                    glob_i =  ((nprocs_per_dim - 1) * N) + i;
                    for (int j = 0; j < N; j++) {
                        glob_j = j;
                        global_output[glob_i][glob_j] = C_n[i][j];
                    }
                }

                // Every other processor's contribution to the global output
                for (int x = 1; x < nprocs; x++) {
                    MPI_Recv(contiguous_local_buffer, N*N, MPI_DOUBLE, x, 0, comm2d, MPI_STATUS_IGNORE);
                    for (int i = 0; i < N; i++) {
                        glob_i =  ((nprocs_per_dim - 1 - floor(x / nprocs_per_dim)) * N) + i;
                        for (int j = 0; j < N; j++) {
                            glob_j = ((x % nprocs_per_dim) * N) + j;
                            global_output[glob_i][glob_j] = contiguous_local_buffer[contiguous_memory_index(i, j, N)];
                        }
                    }
                }

                // Write to file
                write_to_file(global_output, "./final-version/simulation_NTby2_timesteps.txt", nprocs_per_dim, N_glob);
            }
            // For other processors
            else {
                // Initialize contiguous local buffer for each processor's local output to send
                double* contiguous_C_n = contiguous_memory_alloc(N);
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        contiguous_C_n[contiguous_memory_index(i, j, N)] = C_n[i][j];
                    }
                }
                MPI_Send(contiguous_C_n, N*N, MPI_DOUBLE, 0, 0, comm2d);
            }
        } 

        // Write to file upon completion of simulation (time step = NT)
        else if (n == NT - 1) {
            
            // Initialize local buffer for each processor
            double* contiguous_local_buffer = contiguous_memory_alloc(N);

            // For processor 0
            if (mype == 0) {
                // Processor 0's contribution to the global output
                int glob_i;
                int glob_j; 
                for (int i = 0; i < N; i++) {
                    glob_i =  ((nprocs_per_dim - 1) * N) + i;
                    for (int j = 0; j < N; j++) {
                        glob_j = j;
                        global_output[glob_i][glob_j] = C_n[i][j];
                    }
                }

                // Every other processor's contribution to the global output
                for (int x = 1; x < nprocs; x++) {
                    MPI_Recv(contiguous_local_buffer, N*N, MPI_DOUBLE, x, 0, comm2d, MPI_STATUS_IGNORE);
                    for (int i = 0; i < N; i++) {
                        glob_i =  ((nprocs_per_dim - 1 - floor(x / nprocs_per_dim)) * N) + i;
                        for (int j = 0; j < N; j++) {
                            glob_j = ((x % nprocs_per_dim) * N) + j;
                            global_output[glob_i][glob_j] = contiguous_local_buffer[contiguous_memory_index(i, j, N)];
                        }
                    }
                }

                // Write to file
                write_to_file(global_output, "./final-version/simulation_NT_timesteps.txt", nprocs_per_dim, N_glob);
            }
            // For other processors 
            else {
                // Initialize contiguous local buffer for each processor's local output to send
                double* contiguous_C_n = contiguous_memory_alloc(N);
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        contiguous_C_n[contiguous_memory_index(i, j, N)] = C_n[i][j];
                    }
                }
                MPI_Send(contiguous_C_n, N*N, MPI_DOUBLE, 0, 0, comm2d);
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
    int num_threads = stoi(argv[7]);
    string scheme = argv[8];

    cout << "Simulation Parameters:" << endl;
    cout << "N = " << N_glb << endl;
    cout << "NT = " << NT << endl;
    cout << "L = " << L << endl;
    cout << "T = " << T << endl;
    cout << "u = " << u << endl;
    cout << "v = " << v  << endl;
    cout << "Scheme: " << scheme << "\n" << endl;

    // Esimate memory usage
    cout << "Estimated memeory usage = " << N_loc*N_loc*sizeof(double)/1e6 << " MB" << "\n" << endl;

    // Set number of threads
    // num_threads = N_THREADS;
    cout << "Number of threads = " << num_threads << "\n" << endl;
    omp_set_num_threads(num_threads);

    // Start timer
    double t1 = omp_get_wtime();

    // Run simulation
    cout << "Simulating..." << endl;
    advection_simulation(N_loc, NT, L, T, u, v, up, down, left, right, comm2d, mype, nprocs_per_dim, nprocs, scheme);
    cout << "Simulation Complete!" << endl;

    // Stop timer
    double t2 = omp_get_wtime();

    // MPI finalization
    MPI_Finalize();

    // Print time taken
    cout << "Time taken = " << t2-t1 << " seconds" << "\n" << endl;
}