import numpy as np
import matplotlib.pyplot as plt
import ast

def read_file(filename):
    '''
    Read the file and return the data.
    '''
    with open(filename, "r") as f:
        file_str = f.read()
        data = ast.literal_eval(file_str)
        n_procs_per_dim = data[0][0]
        grid = data[1]
        grid = np.array(grid)
    return (n_procs_per_dim, grid)


def generate_plots(filename):
    '''
    Plot and save the data.
    '''
    n_procs_per_dim, data = read_file(filename)
    plt.imshow(data)

    plt.xticks([])
    plt.yticks([])

    # plt.xlim(1, data.shape[0]*n_procs_per_dim)
    # plt.ylim(1, data.shape[0]*n_procs_per_dim)

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar(orientation='vertical')
    plt.title('Advection Simulation Plot')

    plt.savefig('.' + filename.split('.')[1] + '.png')
    
    plt.close()

    


if __name__ == "__main__":
    generate_plots("./final-version/initial_gaussian.txt")

    generate_plots("./final-version/simulation_NTby2_timesteps.txt")

    generate_plots("./final-version/simulation_NT_timesteps.txt")
