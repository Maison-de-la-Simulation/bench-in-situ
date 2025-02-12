# %%
import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="Missing --input argument, please specify input file (txt)")
args=parser.parse_args()

if (args.input):
    toRead = args.input
else:
    print("No input file is given.")
    sys.exit()

def plot_data(data, metrics, cube_size, constraint, pdi_version):
    fig, ax = plt.subplots()
    for item in data:
        if item['metrics'] == metrics:
            # Plot the bar for the current item
            ax.bar(item['label'], item['value'], color='C0')
            ax.set_xlabel('Number of subsets (GCD)')
            ax.set_ylabel("Seconds" if item['unit'] == "s" else item['unit'])
            ax.set_title((f'{metrics} for problem of size {cube_size}, with PDI version {pdi_version}, on {constraint}').replace("_", " "))
    return fig

constraint = ""
pdi_version = ""

def read_data_with_blanks(filename):
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter=' ')
        data = []
        # for row in reader:
        for index, row in enumerate(reader):
            if index == 0:
                constraint = row[-1].split('=')[-1]
                continue
            if index == 1:
                pdi_version = row[0]
                continue
            path, metrics, value, unit = row
            # Label is the simulation size, number of subsets/gcd, and cube size is the problem size
            data.append({'cube_size': path.split('/')[0], 'label': path.split('.')[0].split('s')[1], 'metrics': metrics, 'value': float(value), 'unit': unit})
        # Sort the data based on the integer value of the 'cube_size' and 'label' fields
        data.sort(key=lambda x: (int(x['cube_size']), int(x['label'])))
    return [data, constraint, pdi_version]

# Get unique metrics from the data
metrics = set([item['metrics'] for item in read_data_with_blanks(toRead)[0]])

# Get unique cube sizes from the data
cube_sizes = set([item['cube_size'] for item in read_data_with_blanks(toRead)[0]])
fig, axs = plt.subplots(len(cube_sizes), len(metrics), figsize=(20, 20 * len(cube_sizes)))

# Create a PdfPages object
pdf_pages = PdfPages(datetime.now().strftime("%d_%m_%Y-%H_%M_%S") + ".pdf")

if len(axs.shape) == 1:
    for i, cs in enumerate(cube_sizes):
        for j, m in enumerate(metrics):
            datalist = read_data_with_blanks(toRead)
            fig = plot_data([item for item in datalist[0] if item['cube_size'] == cs], m, cs, datalist[1], datalist[2])
            pdf_pages.savefig(fig)
            axs[j].legend([]) # Clear the legend for cleaner plots
            axs[j].grid(False) # Turn off grid lines for cleaner plots

else:
# Create separate plots for each unique cube size and metric
    for i, cs in enumerate(cube_sizes):
        for j, m in enumerate(metrics):
            datalist = read_data_with_blanks(toRead)
            fig = plot_data([item for item in datalist[0] if item['cube_size'] == cs], m, cs, datalist[1], datalist[2])
            pdf_pages.savefig(fig)
            axs[i, j].legend([]) # Clear the legend for cleaner plots
            axs[i, j].grid(False) # Turn off grid lines for cleaner plots

pdf_pages.close()

# %%
