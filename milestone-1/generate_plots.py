import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
    '''
    Read the file and return the data as a numpy array.
    '''
    data = []
    with open(filename) as f:
        for row in f:
            row_data = []
            for element in row.split(" "):
                row_data.append(np.float64(element))
            data.append(row_data)
    return np.array(data)

def plot_data(filename):
    '''
    Plot and save the data.
    '''
    data = read_file(filename)
    plt.imshow(data)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar(orientation='vertical')
    plt.title('Advection Simulation Plot')
    plt.savefig('.' + filename.split('.')[1] + '.png')
    plt.close()


# Visualize and save data for different timestamps in the simulation
plot_data("./milestone-1/initial_gaussian.txt")

plot_data("./milestone-1/simulation_10000_timesteps.txt")

plot_data("./milestone-1/simulation_15000_timesteps.txt")

plot_data("./milestone-1/simulation_user_specified_timestep.txt")