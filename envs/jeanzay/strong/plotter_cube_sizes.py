# %%
import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="Missing --input argument, please specify input file (txt).")
args=parser.parse_args()

if (args.input):
    toRead = args.input
else:
    print("No input file is given.")
    sys.exit()

# toRead = 'bench_output_MI300_2025-03-26-08-52.txt'

def plot_data(data, metrics, cube_size, constraint, pdi_version):
    fig, ax = plt.subplots()
    for item in data:
        if item['metrics'] == metrics:
            # Plot the bar for the current item
            ax.bar(item['label'], item['value'], color='C0')
            ax.set_xlabel('Number of GPUs')
            ax.set_ylabel("Seconds" if item['unit'] == "s" else item['unit'])
            ax.set_title((f'{metrics} for problem of size {cube_size}, with PDI version {pdi_version}\n'.replace("_", " "))+(f'on {constraint}'))
    return fig

def read_data_with_blanks(filename):
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter=' ')
        data = []
        # for row in reader:
        for index, row in enumerate(reader):
            # Constraint is the name of test case
            if index == 0:
                constraint = row[-1].split('=')[-1]
                continue
            if index == 1:
                pdi_version = row[0]
                continue
            path, metrics, value, unit = row
            if path.startswith("./"):
                path = path[2:]
            # label is the simulation size, number of gpus, and cube size is the problem size
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
pdf_pages = PdfPages("Strong_scaling_plot_" + toRead.split(".txt")[0] + ".pdf")

datalist = read_data_with_blanks(toRead)
constraint = ""
pdi_version = ""

if len(axs.shape) == 1:
# Create separate plots for each metric
    for i, cs in enumerate(cube_sizes):
        for j, m in enumerate(metrics):
                # datalist[1] is the constraint, datalist[2] is the version of pdi
                fig = plot_data([item for item in datalist[0] if item['cube_size'] == cs], m, cs, datalist[1], datalist[2])
                pdf_pages.savefig(fig)
                axs[j].legend([]) # Clear the legend for cleaner plots
                axs[j].grid(False) # Turn off grid lines for cleaner plots

        # Create a new plot for the wall time metrics multiplied by the number of subsets
        wall_time_data = [item for item in datalist[0] if item['cube_size'] == cs and item['metrics'] == 'Wall_time']
        wall_time_values = [item['value'] * int(item['label']) for item in wall_time_data]
        wall_time_labels = [f"{item['label']}" for item in wall_time_data]

        fig, ax = plt.subplots()
        ax.bar(wall_time_labels, wall_time_values, color='C1')
        ax.set_xlabel('Number of GPUs)')
        ax.set_ylabel("Seconds")
        ax.set_title(f'Wall time X Nb GPUs for problem of size {cs}, with PDI version {datalist[2]} \n on {datalist[1]}')
        pdf_pages.savefig(fig)

else:
# Create separate plots for each unique cube size and metric
    for i, cs in enumerate(cube_sizes):
        for j, m in enumerate(metrics):
            fig = plot_data([item for item in datalist[0] if item['cube_size'] == cs], m, cs, datalist[1], datalist[2])
            pdf_pages.savefig(fig)
            axs[i, j].legend([]) # Clear the legend for cleaner plots
            axs[i, j].grid(False) # Turn off grid lines for cleaner plots
        
        # Create a new plot for the wall time metrics multiplied by the number of subsets
        wall_time_data = [item for item in datalist[0] if item['cube_size'] == cs and item['metrics'] == 'Wall_time']
        wall_time_values = [item['value'] * int(item['label']) for item in wall_time_data]
        wall_time_labels = [f"{item['label']}" for item in wall_time_data]

        fig, ax = plt.subplots()
        ax.bar(wall_time_labels, wall_time_values, color='C3')
        ax.set_xlabel('Number of GPUs')
        ax.set_ylabel("Seconds X Number of GPUs")
        ax.set_title(f'Wall time X  Nb GPUs for problem of size {cs}, with PDI version {datalist[2]} \n on {datalist[1]}')
        pdf_pages.savefig(fig)

pdf_pages.close()
plt.close()

# %%
