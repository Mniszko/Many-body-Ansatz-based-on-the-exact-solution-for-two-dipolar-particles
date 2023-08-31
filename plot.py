import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import subprocess
import os

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

def plot_dispersion_output():
    #data = pd.read_csv("dispersion_outputCOM.csv",header=1)
    #plt.plot(data['X'], data['Y'], marker='o',c='orangered',label="center of mass dispersion")

    data = pd.read_csv("dispersion_outputCOG.csv",header=1)
    plt.plot(data['X'], data['Y'], marker='o',c='orchid',label="single particle dispersion")

    plt.xlabel('$L$')
    plt.ylabel(r'$\Delta \rho$')
    
    # Assuming major ticks are at every 10 units (you can adjust this as needed)
    
    major_tick_spacingy = 0.05
    minor_tick_spacingy = major_tick_spacingy / 6
    
    major_tick_spacingx = max(data['X'])/5
    minor_tick_spacingx = major_tick_spacingx / 6
    """
    ax = plt.gca()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(major_tick_spacingx))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(minor_tick_spacingx))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(major_tick_spacingy))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(minor_tick_spacingy))
    ax.grid(which='minor', linestyle='--', linewidth=0.5)
    """
    plt.legend()
    plt.title('Dispersion of lowest energy eigenfunction for $g=1$, Dirac delta interaction')
    plt.grid(True, which='both')  # This will enable grid for both major and minor ticks
    plt.show()

def plot_g_output():
    data = pd.read_csv("g_output.csv",header=1)
    
    plt.plot(data['X'], data['Y'], marker='o',c='orangered')
    plt.xlabel('$g$')
    plt.ylabel('$E_0$')
    plt.title('$E_0$ for $L=1$, Dirac delta interaction')
    plt.grid(True)
    plt.show()

def plot_first_coefficient():
    directory="./vector_outputs/"
    # 1. List all files in the directory that match the pattern "lowest_wavefunction*.csv".
    all_files = [f for f in os.listdir(directory) if f.startswith("lowest_wavefunction") and f.endswith(".csv")]
    
    # 2. Extract the floating point number <n> from each filename.
    def extract_number(filename):
        # Extract the number between "lowest_wavefunction" and ".csv"
        start_index = len("lowest_wavefunction")
        end_index = filename.find(".csv")
        return float(filename[start_index:end_index])
    
    # 3. Sort the filenames based on the extracted numbers from lowest to highest.
    all_files.sort(key=extract_number)
    
    n_values = []
    data_values = []
    
    # 4. For each file, open it, read the last column, and extract the first element.
    for file in all_files:
        filepath = os.path.join(directory, file)
        df = pd.read_csv(filepath)
        first_element = df.iloc[0, -1]
        
        n_values.append(extract_number(file))
        data_values.append(abs(first_element))
    
    # 5. Plot the numbers against their corresponding <n> values.
    plt.plot(n_values, data_values, '-o')
    plt.xlabel('$L$')
    plt.ylabel('$c_{0,0}$')
    plt.title('plot of first coefficients for lowest energy level state in function of length')
    plt.grid(True)
    plt.show()


if __name__ == "__main__": 
    plot_dispersion_output()
    #plot_length_output()
    #plot_first_coefficient()
    #plot_g_output()
