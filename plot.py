import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter
import matplotlib
from matplotlib.ticker import ScalarFormatter, MultipleLocator, AutoMinorLocator
import subprocess
import locale
import glob
import os
import re
import numpy as np

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

def plot_length_output(loglog=True):    
    
    filenames = glob.glob("length_output[0-9][0-9][0-9].csv")
    print(filenames)
    for fname in filenames:
        data = pd.read_csv(fname, header=1)
        
        numbers = re.findall(r'(\d)', fname)
        legend_label = '({})'.format(', '.join(numbers))
        if loglog:
            plt.loglog(data['X'], data['Y'], marker='o', label=legend_label)  
        else:
            plt.plot(data['X'], data['Y'], marker='o', label=legend_label)
    

    #reference line - 1D fermionization energy
    refline = lambda L : 2+4*np.pi**2/L**2*(2-1/2)/6
    #refline = lambda L : 2+4*np.pi**2/L**2
    L = np.linspace(1e-6,20,int(2e3))
    plt.plot(L,refline(L),ls="dashed",label="analytic solution in 1D")



    plt.ylim(2,5)
    plt.xlim(-1e-1,20)
    plt.xlabel('$L$')
    plt.ylabel('$E_0$')
    plt.legend()
    #plt.title('$E_0$ for $g=1$, Dirac delta interaction')
    plt.grid(True)
    plt.show()

def plot_dispersion_output():
    data = pd.read_csv("dispersion_outputCOM.csv",header=1)
    plt.plot(data['X'], data['Y'], marker='o',c='orangered',label="center of mass dispersion")

    data = pd.read_csv("dispersion_outputCOG.csv",header=1)
    plt.plot(data['X'], data['Y'], marker='o',c='orchid',label="single particle dispersion")

    plt.xlabel('$L$')
    plt.ylabel(r'$\Delta \rho$')
    
    # Assuming major ticks are at every 10 units (you can adjust this as needed)
    
    major_tick_spacingy = 0.05
    minor_tick_spacingy = major_tick_spacingy / 6
    
    major_tick_spacingx = max(data['X'])/5
    minor_tick_spacingx = major_tick_spacingx / 6
    
    ax = plt.gca()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(major_tick_spacingx))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(minor_tick_spacingx))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(major_tick_spacingy))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(minor_tick_spacingy))
    ax.grid(which='minor', linestyle='--', linewidth=0.5)
    
    plt.xlim(-5e-1,8)
    plt.plot(np.array([-1,8]),[0.42625,0.42625],c='brown',ls="dashed")
    plt.plot(np.array([-1,8]),[0.30140,0.30140],c='brown',ls="dashed",label=r"dispersion in $L\rightarrow\infty$")
    plt.legend()
    #plt.title('Dispersion of lowest energy eigenfunction for $g=1$, Dirac delta interaction')
    plt.grid(True, which='both')  # This will enable grid for both major and minor ticks
    plt.show()

def plot_g_output():
    
    filenames = glob.glob("g_output[0-9][0-9][0-9].csv")
    print(filenames)
    for fname in filenames:
        data = pd.read_csv(fname, header=1)
        
        numbers = re.findall(r'(\d)', fname)
        legend_label = '({})'.format(', '.join(numbers))
        
        plt.plot(data['X'], data['Y'], marker='o', label=legend_label) 
    
    xx = np.linspace(0,20,20)
    yy = np.zeros(20)
    #yy += np.pi**2/3
    #yy += 2+4*np.pi**2/1*(2-1/2)/6
    #plt.plot(xx,yy)
    
    plt.xlabel('$g$')
    plt.ylabel('$E_0$')
    #plt.title('$E_0$ for $L=1$, Dirac delta interaction')
    #plt.legend()
    plt.grid(True)
    plt.show()
    plt.savefig("delta-g-plot.svg",format="svg")

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
    #plt.title('plot of first coefficients for lowest energy level state in function of length')
    plt.grid(True)
    plt.show()

if __name__ == "__main__": 
    scale = 6
    plt.rc('text', usetex=False)
    plt.rcParams['figure.figsize']=[8,6]
    plt.rc('font', size=10*scale/2)          # controls default text sizes
    plt.rc('axes', titlesize=10*scale/3.5)     # fontsize of the axes title
    plt.rc('axes', labelsize=10*scale/3.5)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=8*scale/3.5)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=8*scale/3.5)    # fontsize of the tick labels
    plt.rc('legend', fontsize=8*scale/4)    # legend fontsize
    plt.rc('figure', titlesize=16*scale/3.5)
    locale.setlocale(locale.LC_NUMERIC, "pl_PL.UTF-8")
    plt.rcParams['axes.formatter.use_locale']=True
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'


    # FUNKCJE NIEAKTUALNE! nie uwzględniają nowych nazw plików
    plot_length_output()
    plot_length_output(False)
    #plot_first_coefficient()
    
    #plot_g_output()
    #plot_dispersion_output()
