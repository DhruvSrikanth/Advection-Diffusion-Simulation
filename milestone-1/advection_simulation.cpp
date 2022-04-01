#include <iostream>
using namespace std;
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>

vector<vector<double> > advection_simulation(int N, int NT, double L, double T, double u, double v) {
    
    // Initialize variables
    vector<vector<double> > C_n(N, vector<double>(N, 0.0));
    vector<vector<double> > C_n_1(N, vector<double>(N, 0.0));

    double dx = L/N;
    double dt = T/NT;
    
    // Check Courant stability condition
    assert(dt<=dx/sqrt(2*(u*u + v*v)));

    // Initialize gaussian grid
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C_n[i][j] = exp(-((i-L/2)*(i-L/2) + (j-L/2)*(j-L/2))/(2*L/16));
        }
    }


    // Compute next time step
    for (int n = 0; n < NT; n ++) {
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
        // Update C_n
        for (int i = 0; i < N; i ++) {
            for (int j = 0; j < N; j ++) {
                C_n[i][j] = C_n_1[i][j];
            }
        }
    }
    return C_n;
}

void write_to_file(vector<vector<double> > C_n, int N, string filename) {
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
    vector<vector<double> > C_n = advection_simulation(N, NT, L, T, u, v);
    cout << "Simulation Complete!" << endl;

    // Write output to file
    cout << "Writing output to file..." << endl;
    write_to_file(C_n, N, "./milestone-1/simulation_user_specified_timestep.txt");
    cout << "Output written to file!" << "\n" << endl; 

    // Get initial gaussian
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C_n[i][j] = exp(-((i-L/2)*(i-L/2) + (j-L/2)*(j-L/2))/(2*L/16));
        }
    }

    // Write output to file
    cout << "Writing initial gaussian to file..." << endl;
    write_to_file(C_n, N, "./milestone-1/initial_gaussian.txt");
    cout << "Initial gaussian written to file!" << "\n" << endl;

    // Perform simulation
    NT = 10000;
    cout << "Simulating for 10000 timesteps..." << endl;
    C_n = advection_simulation(N, NT, L, T, u, v);
    cout << "Simulation Complete!" << endl;

    // Write output to file
    cout << "Writing output to file..." << endl;
    write_to_file(C_n, N, "./milestone-1/simulation_10000_timesteps.txt");
    cout << "Output written to file!" << "\n" << endl;

    // Perform simulation
    NT = 15000;
    cout << "Simulating for 15000 timesteps..." << endl;
    C_n = advection_simulation(N, NT, L, T, u, v);
    cout << "Simulation Complete!" << endl;

    // Write output to file
    cout << "Writing output to file..." << endl;
    write_to_file(C_n, N, "./milestone-1/simulation_15000_timesteps.txt");
    cout << "Output written to file!" << "\n" << endl;


    

    return 0;

}