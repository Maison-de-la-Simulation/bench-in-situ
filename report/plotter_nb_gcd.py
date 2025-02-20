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

# toRead = 'ALL-output_v4_noGenoa.txt'

def plot_data(data, metrics, label, constraint, pdi_version):
    fig, ax = plt.subplots()
    for item in data:
        if item['metrics'] == metrics:
            # Plot the bar for the current item
            ax.bar(item["cube_size"], item['value'], color='C0')
            ax.set_xlabel('Problem size')
            ax.set_ylabel("Seconds" if item['unit'] == "s" else item['unit'])
            ax.set_title((f'{metrics} with {item["label"]} GCD, with PDI version {pdi_version}, on {constraint}').replace("_", " "))
    return fig

def read_data_with_blanks(filename):
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter=' ')
        data = []
        # for row in reader:
        for index, row in enumerate(reader):
            # Constraint is the GCD type
            if index == 0:
                constraint = row[-1].split('=')[-1]
                continue
            if index == 1:
                pdi_version = row[0]
                continue
            path, metrics, value, unit = row
            if path.startswith("./"):
                path = path[2:]
            # label is the simulation size, number of subsets/gcd, and cube size is the problem size
            data.append({'cube_size': path.split('/')[0], 'label': path.split('.')[0].split('s')[1], 'metrics': metrics, 'value': float(value), 'unit': unit})
        # Sort the data based on the integer value of the 'cube_size' and 'label' fields
        data.sort(key=lambda x: (int(x['cube_size']), int(x['label'])))
    return [data, constraint, pdi_version]

# Get unique metrics from the data
metrics = set([item['metrics'] for item in read_data_with_blanks(toRead)[0]])

# Get unique cube sizes from the data
cube_sizes = set([item['cube_size'] for item in read_data_with_blanks(toRead)[0]])

nb_GCD_sizes = set([item['label'] for item in read_data_with_blanks(toRead)[0]])

# Create a PdfPages object
pdf_pages = PdfPages("Number of GCD plot - " + toRead + ".pdf")

datalist = read_data_with_blanks(toRead)
constraint = ""
pdi_version = ""

fig, axs = plt.subplots(len(nb_GCD_sizes), len(metrics), figsize=(20, 20 * len(nb_GCD_sizes)))

if len(axs.shape) == 1:
# Create separate plots for each metric
    for i, lb in enumerate(nb_GCD_sizes):
        for j, m in enumerate(metrics):
                # datalist[1] is the constraint, datalist[2] is the version of pdi
                fig = plot_data([item for item in datalist[0] if item['label'] == lb], m, lb, datalist[1], datalist[2])
                pdf_pages.savefig(fig)
                axs[j].legend([]) # Clear the legend for cleaner plots
                axs[j].grid(False) # Turn off grid lines for cleaner plots

        # # Create a new plot for the wall time metrics multiplied by the number of subsets
        # wall_time_data = [item for item in datalist[0] if item['label'] == lb and item['metrics'] == 'Wall_time']
        # wall_time_values = [item['value'] * int(item['cube_size']) for item in wall_time_data]
        # wall_time_labels = [f"{item['cube_size']}" for item in wall_time_data]

        # fig, ax = plt.subplots()
        # ax.bar(wall_time_labels, wall_time_values, color='C3')
        # ax.set_xlabel('Problem size')
        # ax.set_ylabel("Seconds X Problem size")
        # ax.set_title(f'Wall time X problem size for number of GCD {lb}, with PDI version {datalist[2]}, on {datalist[1]}')
        # pdf_pages.savefig(fig)

else:
# Create separate plots for each unique cube size and metric
    for i, lb in enumerate(nb_GCD_sizes):
        for j, m in enumerate(metrics):
            fig = plot_data([item for item in datalist[0] if item['label'] == lb], m, lb, datalist[1], datalist[2])
            pdf_pages.savefig(fig)
            axs[i, j].legend([]) # Clear the legend for cleaner plots
            axs[i, j].grid(False) # Turn off grid lines for cleaner plots
        
        # # Create a new plot for the wall time metrics multiplied by the number of subsets
        # wall_time_data = [item for item in datalist[0] if item['label'] == lb and item['metrics'] == 'Wall_time']
        # wall_time_values = [item['value'] * int(item['cube_size']) for item in wall_time_data]
        # wall_time_labels = [f"{item['cube_size']}" for item in wall_time_data]

        # fig, ax = plt.subplots()
        # ax.bar(wall_time_labels, wall_time_values, color='C3')
        # ax.set_xlabel('Problem size')
        # ax.set_ylabel("Seconds X Problem size")
        # ax.set_title(f'Wall time X problem size for number of GCD {lb}, with PDI version {datalist[2]}, on {datalist[1]}')
        # pdf_pages.savefig(fig)

pdf_pages.close()
plt.close()

# %%
