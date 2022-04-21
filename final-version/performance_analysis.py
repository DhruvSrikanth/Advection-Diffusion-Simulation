import numpy as np
import matplotlib.pyplot as plt


def generate_strong_scaling_plot(n, time, analysis_num):
    """
    Generates a plot of the strong scaling results.
    """
    t1 = time[0]
    Sn = [t1/t for t in time]

    n = np.array(n)
    Sn = np.array(Sn)

    # Plot the data
    plt.plot(n, Sn, '-o')
    # Set the labels
    plt.xlabel("n (number of nodes)")
    plt.ylabel("S(n)")
    # Set the title
    plt.title("Strong Scaling Analysis")
    # Save the plot
    plt.savefig("./final-version/strong_scaling_{}.png".format(analysis_num))
    plt.clf()

def generate_weak_scaling_plot(n, time, analysis_num):
    """
    Generates a plot of the weak scaling results.
    """
    t1 = time[0]
    Sn = [t1/t for t in time]

    n = np.array(n)
    Sn = np.array(Sn)

    # Plot the data
    plt.plot(n, Sn, '-o')
    # Set the labels
    plt.xlabel("n (number of nodes)")
    plt.ylabel("S(n)")
    # Set the title
    plt.title("Weak Scaling Analysis")
    # Save the plot
    plt.savefig("./final-version/weak_scaling_{}.png".format(analysis_num))
    plt.clf()

n1 = [1,4,16]
time1 = [416.7, 77.8, 11.6]
generate_strong_scaling_plot(n1, time1, "1")

n2 = n1
time2 = [232.6, 47.5, 11.7]
generate_strong_scaling_plot(n2, time2, "2")

n3 = n1
time3 = [384.6, 400.0, 425.5]
generate_weak_scaling_plot(n3, time3, "1")


n4 = n1
time3 = [434.8, 434.8, 454.5]
generate_weak_scaling_plot(n3, time3, "2")
