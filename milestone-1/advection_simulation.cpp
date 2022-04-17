#include <iostream>
using namespace std;
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <chrono>
using namespace std::chrono;

void write_to_file(double** const& C_n, int const& N, string const& filename) {
    ofstream out(filename);
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            if (j == N-1) {
                out << C_n[i][j] << "\n";
            } else {
                out << C_n[i][j] << " ";
            }
        }
    }
    out.close();
}

double** create_matrix(int const& N){
    double** mat = new double*[N];

    for (int i = 0; i < N; i++) {
        mat[i] = new double[N];
        for (int j = 0; j < N; j++) {
            mat[i][j] = 0.0;
        }
    }
    return mat;
}

void initial_gaussian(double**& C_n, int const& N, double const& L) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C_n[i][j] = exp(-(pow((i-N/2)*L/N, 2) + pow((j-N/2)*L/N, 2))/(L*L/8));
        }
    }
}

void apply_boundary_conditions(int const& N, double const& dt, double const& dx, double const& u, double const& v, double** const& C_n, double**& C_n_1) {
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            double C_n_up;
            double C_n_down;
            double C_n_left;
            double C_n_right;

            // Apply periodic boundary conditions
            if (i == 0) {
                C_n_up = C_n[N-1][j];
                C_n_down = C_n[i+1][j];
                if (j == 0) {
                    C_n_left = C_n[i][N-1];
                    C_n_right = C_n[i][j+1];
                } else if (j == N-1) {
                    C_n_left = C_n[i][j-1];
                    C_n_right = C_n[i][0];
                } else {
                    C_n_left = C_n[i][j-1];
                    C_n_right = C_n[i][j+1];
                }
            } 
            else if (i == N-1) {
                C_n_up = C_n[i-1][j];
                C_n_down = C_n[0][j];
                if (j == 0) {
                    C_n_left = C_n[i][N-1];
                    C_n_right = C_n[i][j+1];
                } else if (j == N-1) {
                    C_n_left = C_n[i][j-1];
                    C_n_right = C_n[i][0];
                } else {
                    C_n_left = C_n[i][j-1];
                    C_n_right = C_n[i][j+1];
                }
            }
            else {
                C_n_up = C_n[i-1][j];
                C_n_down = C_n[i+1][j];
                if (j == 0) {
                    C_n_left = C_n[i][N-1];
                    C_n_right = C_n[i][j+1];
                } else if (j == N-1) {
                    C_n_left = C_n[i][j-1];
                    C_n_right = C_n[i][0];
                } else {
                    C_n_left = C_n[i][j-1];
                    C_n_right = C_n[i][j+1];
                } 
            }
            
            // Compute next time step
            C_n_1[i][j] = 0.25*(C_n_up + C_n_down + C_n_left + C_n_right) - (dt/(2*dx))*(u*(C_n_down - C_n_up) + v*(C_n_right - C_n_left));
        }
    }
}

void advection_simulation(int N, int NT, double L, double T, double u, double v) {
    
    // Initialize variables
    double** C_n = create_matrix(N);
    double** C_n_1 = create_matrix(N);

    double dx = L/N;
    double dt = T/NT;

    // Check Courant stability condition
    assert(dt<=dx/sqrt(2*(u*u + v*v)));

    // Initialize gaussian grid
    initial_gaussian(C_n, N, L);
    // Apply boundary conditions
    apply_boundary_conditions(N, dt, dx, u, v, C_n, C_n_1);
    // Update C_n
    swap(C_n, C_n_1);

   // Write output to file
    cout << "Writing initial gaussian to file..." << endl;
    write_to_file(C_n, N, "./milestone-1/initial_gaussian.txt");
    cout << "Initial gaussian written to file!" << "\n" << endl;

    // Compute next time step
    for (int n = 0; n < NT; n ++) {
        auto start = high_resolution_clock::now();
        apply_boundary_conditions(N, dt, dx, u, v, C_n, C_n_1);
        swap(C_n, C_n_1);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Iteration: " << n << " - Grind Rate: " << 1e6/duration.count() << " iter/sec" << endl;

        if (n == 10000-1) {
            cout << "Writing output to file..." << endl;
            write_to_file(C_n, N, "./milestone-1/simulation_10000_timesteps.txt");
            cout << "Output written to file!" << "\n" << endl; 
        }
        else if (n == 15000-1) {
            cout << "Writing output to file..." << endl;
            write_to_file(C_n, N, "./milestone-1/simulation_15000_timesteps.txt");
            cout << "Output written to file!" << "\n" << endl; 
        }
    }

    // Write output to file
    cout << "Writing output to file..." << endl;
    write_to_file(C_n, N, "./milestone-1/simulation_user_specified_timestep.txt");
    cout << "Output written to file!" << "\n" << endl; 
}

int main(int argc, char** argv) {
    // Initialize variables
    int N = stoi(argv[1]);
    int NT = stoi(argv[2]);
    double L = stod(argv[3]);
    double T = stod(argv[4]);
    double u = stod(argv[5]);
    double v = stod(argv[6]);

    cout << "Simulation Parameters:" << endl;
    cout << "N = " << N << endl;
    cout << "NT = " << NT << endl;
    cout << "L = " << L << endl;
    cout << "T = " << T << endl;
    cout << "u = " << u << endl;
    cout << "v = " << v << "\n" << endl;

    // Esimate memory usage
    cout << "Estimated memeory usage = " << N*N*sizeof(double)/1e6 << " MB" << "\n" << endl;

    // Perform simulation
    cout << "Simulating..." << endl;
    auto s = high_resolution_clock::now();
    advection_simulation(N, NT, L, T, u, v);
    auto e = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(e - s);
    cout << 1e6*duration.count() << " sec" << endl;
    cout << "Simulation Complete!" << endl;

    return 0;

}