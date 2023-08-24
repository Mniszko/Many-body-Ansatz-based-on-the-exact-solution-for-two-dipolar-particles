import pandas as pd
import matplotlib.pyplot as plt
import subprocess

def plot_csv_data(filename):
    # Read the CSV file using pandas
    subprocess.run(["make"])
    subprocess.run(["./g-plot"])

    data = pd.read_csv(filename)
    
    # Plot the data
    plt.plot(data['X'], data['Y'], marker='o')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Plot of Y vs X')
    plt.grid(True)
    plt.show()

def plot_length_output():
    data = pd.read_csv("length_output.csv",header=1)
    
    plt.plot(data['X'], data['Y'], marker='o',c='orangered')
    plt.xlabel('$L$')
    plt.ylabel('$E_0$')
    plt.title('$E_0$ for $g=1$, Dirac delta interaction')
    plt.grid(True)
    plt.show()

def plot_g_output():
    data = pd.read_csv("g_output.csv",header=1)
    
    plt.plot(data['X'], data['Y'], marker='o',c='orangered')
    plt.xlabel('$g$')
    plt.ylabel('$E_0$')
    plt.title('$E_0$ for $L=1$, Dirac delta interaction')
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    plot_length_output()
    plot_g_output()
