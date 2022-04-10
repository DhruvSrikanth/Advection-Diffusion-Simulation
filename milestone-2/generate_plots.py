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
    plt.xlabel("n (number of cores)")
    plt.ylabel("S(n)")
    # Set the title
    plt.title("Strong Scaling Analysis")
    # Save the plot
    plt.savefig("./milestone-2/strong_scaling_{}.png".format(analysis_num))
    plt.clf()

def generate_weak_scaling_plot(n, time):
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
    plt.xlabel("n (number of cores)")
    plt.ylabel("S(n)")
    # Set the title
    plt.title("Weak Scaling Analysis")
    # Save the plot
    plt.savefig("./milestone-2/weak_scaling.png")
    plt.clf()

n1 = [1,2,3,4,5,6,7,8]
time1 = [36.4707, 18.6834, 12.7379, 9.7075, 7.7995, 7.3461, 7.1528, 7.5088]
generate_strong_scaling_plot(n1, time1, "1")

n2 = n1
time2 = [0.15496, 0.096304, 0.090706, 0.058696, 0.064069, 0.068613, 0.087449, 0.092364]
generate_strong_scaling_plot(n2, time2, "2")

n3 = n1
time3 = [2.34124, 2.35585, 2.41521, 2.40799, 2.4846, 2.62962, 3.41103, 3.96412]
generate_weak_scaling_plot(n3, time3)


