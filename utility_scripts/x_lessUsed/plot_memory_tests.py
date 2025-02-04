#!/usr/bin/env python3

# run this in directory with memory tracking text files

import os
import glob
import re
import matplotlib.pyplot as plt

output_dir = os.getcwd()
os.chdir(output_dir)

postfix_list = ['K','M','G']

mem_list = []
run_name_list = []

for filename in glob.glob("run_memory_info_*.txt"):
    
    try:
        found = re.search('run_memory_info_(.+?)\.txt', filename).group(1)
    except AttributeError:
        found = "missing_run_metadata"

    run_name_list.append(found)
    with open(filename) as file:
        lines = [line.rstrip() for line in file]
        good_lines = []
        for line in lines:
            if len(line) > 0:
                postfix = line[-1]
                if postfix in postfix_list:
                    value = float(line[:-1])
                    if postfix == 'K':
                        good_lines.append(value/1000000)
                    elif postfix == 'M':
                        good_lines.append(value/1000)
                    else:
                        good_lines.append(value)
        mem_list.append(good_lines)

fig,ax = plt.subplots()

for ii in range(len(mem_list)):
    ax.plot(mem_list[ii],label = run_name_list[ii])
    
ax.legend()
ax.set_xlabel("Runtime (x 10 seconds)")
ax.set_ylabel("Memory usage (GB)")

plt.title('Memory usage (GB) vs Runtime (not simulation time)')

plt.show()

