#!/usr/bin/env python3

# run this in directory with log files

import os
import glob
import re
import matplotlib.pyplot as plt

output_dir = os.getcwd()
os.chdir(output_dir)

runtime_string = 'total runtime:'
exectime_string = 'total execution time:'

run_name_list = []
file_name_list = []
calc_length_list = []
runtime_list = []
exectime_list = []

for filename in glob.glob("log_*.txt"):
    
    try:
        found = re.search('log_(.+?)\.txt', filename).group(1)
    except AttributeError:
        found = "missing_run_metadata"

    file_name_list.append(filename)
    run_name_list.append(found)

run_name_list.sort()
file_name_list.sort()

for filename in file_name_list:

    calc_length_pre1 = filename.split('.txt')[0]
    calc_length = int(calc_length_pre1.split('_')[-3])

    calc_length_list.append(calc_length)

    print(calc_length)

    with open(filename) as file:
        for line in file:
            if re.search(runtime_string, line):
                #runtime_list.append(int(float(line.split(':')[1].rstrip())/60))
                runtime_list.append(float(line.split(':')[1].rstrip())/60)
                print(float(line.split(':')[1].rstrip())/60)
            if re.search(exectime_string, line):
                #exectime_list.append(int(float(line.split(':')[1].rstrip())/60))
                exectime_list.append(float(line.split(':')[1].rstrip())/60)

slope = (runtime_list[-1]-runtime_list[0])/(calc_length_list[-1]-calc_length_list[0])

fig, ax = plt.subplots()
ax.plot(calc_length_list[:],runtime_list[:],marker='*', label='run time')
ax.plot(calc_length_list[:],exectime_list[:],marker='*', label='execution time')

ax.set_aspect('equal')

plt.title("Comparing run time and execution time\n for different calc timestep lengths.\n Slope ~= {:.02f}".format(slope))
ax.set_xlabel("Number of seedings")
ax.set_ylabel("Time (minutes)")

plt.legend()

for i_x, i_y in zip(calc_length_list,runtime_list):
        plt.text(i_x, i_y, '({}, {:.01f})'.format(i_x, i_y))

plt.grid()
plt.show()

