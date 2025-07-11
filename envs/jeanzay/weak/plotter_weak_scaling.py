# %%
import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import argparse

def read_data_with_blanks(toRead):
    with open(toRead, 'r') as file:
        reader = csv.reader(file, delimiter=' ')
        data = []
        for index, row in enumerate(reader):
            if index == 0:
                constraint = row[-1].split('=')[-1]
                continue
            if index == 1:
                pdi_version = row[0]
                continue
            path, metrics, value, unit = row
            if path.startswith("./"):
                path = path[2:]
            data.append({'cube_size': path.split('/')[0], 'label': path.split('.')[0].split('s')[1], 'metrics': metrics, 'value': float(value), 'unit': unit})
        data.sort(key=lambda x: (int(x['cube_size']), int(x['label'])))
        return [data, constraint, pdi_version]

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="Missing --input argument, please specify input file (txt).")
args = parser.parse_args()

if args.input:
    toRead = args.input
else:
    print("No input file is given.")
    sys.exit()

# toRead = 'bench_weakScaling_MI300_2025-03-25-09-51.txt'

# Read the data
data, constraint, pdi_version = read_data_with_blanks(toRead)

# Extract the relevant parts of the data
parsed_data = []
for entry in data:
    cube_size = entry['cube_size']
    label = entry['label']
    metrics = entry['metrics']
    value = entry['value']
    unit = entry['unit']
    parsed_data.append({'cube_size': cube_size, 'label': label, 'metric': metrics, 'value': value, 'unit': unit})

# Create a PDF file to save the plots
pdf_pages = PdfPages("Weak_scaling_plot_" + toRead.split(".txt")[0] + ".pdf")

# Function to create a subplot
def create_subplot(metric, data):
    plt.figure(figsize=(10, 6))
    x_values = [entry['label'] for entry in data if entry['metric'] == metric]
    y_values = [entry['value'] for entry in data if entry['metric'] == metric]
    plt.plot(x_values, y_values, label=metric.replace("_", " "))
    plt.xlabel('Nb GPUs')

    ## get units ==> all the units must be the same
    all_units = [entry['unit'] for entry in data if entry['metric'] == metric]
    plt.ylabel("Seconds" if all_units[0] == "s" else all_units[0])

    plt.title((f'{metric}, with PDI version {pdi_version}'.replace("_", " "))+(f', on {constraint}'))
    plt.legend()
    pdf_pages.savefig()  # Save the current figure into the PDF file
    plt.close()  # Close the figure to free up memory

# Performance plot
create_subplot('Performance', parsed_data)

# Wall_time plot
create_subplot('Wall_time', parsed_data)

# I/O_time plot
create_subplot('I/O_time', parsed_data)

# Compute_time plot
create_subplot('Compute_time', parsed_data)

# Close the PDF file
pdf_pages.close()
