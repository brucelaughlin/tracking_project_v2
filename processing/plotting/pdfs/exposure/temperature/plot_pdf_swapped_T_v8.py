# v8: just further updating... variable names, then colors...

# v7: now using .npz files, and np.save (pickle is old and bug prone)

# v6: use SO advice about colorbar

# V5: V3 I actually liked the 2-label color bars on V3, so copy that, and change the y-ticks of the main plot, and maybe add grid lines

# V2: updates...

#------------------------------------------------------------------
#pdf_file_name_pre = "pdf_data_output_seasonal_v3_tenFileTest.p"
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_run_test3.p"
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_test4_physics_only_AKs_1en5.p"
#pdf_file_name_pre = "pdf_data_output_seasonal_ranges_O2_pH___pld_20_30_z_two_file_test.npz"
pdf_file_name_pre = "pdf_data_output_seasonal_ranges_O2_pH___pld_20_29_z_two_file_test_swapped.npz"
#------------------------------------------------------------------

#------------------------------------------------------------------
#pdf_file_name = pdf_file_name_pre[0:-2] + ".p"
pdf_file_name = pdf_file_name_pre
#------------------------------------------------------------------


#------------------------------------------------------------------
fig_paramTitle = "wc15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-day PLD"
fig_mainTitle = "PDFs of average temperature experienced (y-axis) vs release box number (x-axis).\nGrouped according to season of release."
#fig_mainTitle = "PDFs of Average Temperature Experienced (y-axis) vs Release Box Number (x-axis).\nGrouped according to season of release."
#fig_mainTitle = "PDFs of Average Temperature Experienced (y-axis) vs Release Location (x-axis).\nGrouped according to season of release."

fig_fullTitle = fig_mainTitle + "\n" + fig_paramTitle

yLabel = "average temperature experienced (C$^\circ$)"


#------------------------------------------------------------------


import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import datetime

#------------------------------------------------------------------
base_path = '/home/blaughli/tracking_project/'
pdf_raw_directory = base_path + 'practice/bounding_boxes/final_locations/z_output/'

pdf_raw_file = pdf_raw_directory + pdf_file_name
#------------------------------------------------------------------


# Trim histograms... done by eye by plotting this same plot...
T_min = 5
T_max = 25
#n_ticks = int((T_max-T_min)/2) - 1 #kinda a hack, seems to match the naive plot without specified ticks
n_ticks = int((T_max-T_min)/2) + 1  # Now want ticks every 2 degres

bin_size = 0.1
bin_T_min = int(T_min/bin_size)
bin_T_max = int(T_max/bin_size)
bin_T_range = bin_T_max - bin_T_min

y_ticks = np.linspace(0,bin_T_range,n_ticks)
y_ticks = [int(x) for x in y_ticks]

y_tick_labels =  np.linspace(T_min,T_max,n_ticks)
y_tick_labels = [str(x) for x in y_tick_labels]



d = np.load(pdf_raw_file)

pdf_arrays_T = d['pdf_arrays_T']
#pdf_arrays_O2 = d['pdf_arrays_O2']
#pdf_arrays_pH = d['pdf_arrays_pH']
#pdf_arrays_connectivity = d['pdf_arrays_connectivity']
#pdf_arrays_settleTime = d['pdf_arrays_settleTime']
#counter_array = d['counter_array']
#oxygen_limit_list = d['oxygen_limit_list']
#pH_limit_list = d['pH_limit_list']
box_num_mod = d['box_num_mod']
tick_positions = d['tick_positions']
tick_labels = d['tick_labels']
first_continent_box_dex = d['first_continent_box_num']



pdf_list_full = []

pdf_max_val = -999999
pdf_min_val = 999999

# -------------------------------------------------------------------
for pdf in pdf_arrays_T[1:]:
#for pdf in pdf_arrays_T[1:]:
    pdf_full = pdf.copy()
    pdf_full = pdf_full[:,bin_T_min:bin_T_max+1]
    row_sums = pdf_full.sum(axis=1)
    pdf_full = pdf_full / row_sums[:, np.newaxis]
    pdf_full = np.log10(pdf_full)
    pdf_list_full.append(pdf_full)
    if np.amax(pdf_full) > pdf_max_val:
        pdf_max_val = np.amax(pdf_full)
    if np.amin(np.ma.masked_invalid(pdf_full)) < pdf_min_val:
        pdf_min_val = np.amin(np.ma.masked_invalid(pdf_full))
# -------------------------------------------------------------------


num_dummy_lines = 1

n_rows = int(np.shape(pdf_list_full[0])[0]) + num_dummy_lines
n_columns = int(np.shape(pdf_list_full[0])[1])
X = np.arange(-0.5, n_rows, 1)
Y = np.arange(-0.5, n_columns, 1)


#Modify "tick_positions" for plotting (we add an empty bin before the first continent box)
tick_positions_modified = tick_positions.copy()
for ii in range(len(tick_positions_modified)):
    if tick_positions_modified[ii] >= first_continent_box_dex:
        tick_positions_modified[ii] += num_dummy_lines


stagger_dex = 0
tick_labels_double_X = []
for ii in range(len(tick_labels)):
    stagger_dex += 1
    if (tick_positions[ii] < first_continent_box_dex) and (stagger_dex % 2 == 0):
        tick_labels_double_X.append("{}\n\n{}".format(tick_positions[ii]+1,tick_labels[ii]))
    else:
        tick_labels_double_X.append("{}\n{}".format(tick_positions[ii]+1,tick_labels[ii]))




seasons = ["(Dec, Jan, Feb)","(Mar, Apr, May)","(Jun, Jul, Aug)","(Sep, Oct, Nov)"]
season_titles = []
for ii in range(len(seasons)):
    season_titles.append("{}".format(seasons[ii]))


x_line_positions = []
for ii in range(np.shape(pdf_list_full[0])[0]):
    if ii < first_continent_box_dex:
        if ii%2 == 0:
            x_line_positions.append(ii+1.5)
    else:
        x_line_positions.append(ii+0.5)



v_scale = 6

fig,axs = plt.subplots(2,2)


plt.setp(axs,xticks=tick_positions_modified,xticklabels=tick_labels_double_X)


label_fontSize = 10
#label_fontSize = 8

for ii in range(len(pdf_list_full)):

    # 2D full plot (minus day 1)
    pdf_plot = pdf_list_full[ii]
    pdf_separated = np.empty((np.shape(pdf_plot)[0] + num_dummy_lines,np.shape(pdf_plot)[1]))
    pdf_separated[:] = np.nan
    pdf_separated[0:first_continent_box_dex,:] = pdf_plot[0:first_continent_box_dex,:]
    pdf_separated[first_continent_box_dex + num_dummy_lines:,:] = pdf_plot[first_continent_box_dex:,:]



    if ii == 0:
    #if ii == 1:
        mesh1 = axs[0,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,0].title.set_text(season_titles[ii])
        axs[0,0].set_yticks(y_ticks)
        axs[0,0].set_yticklabels(y_tick_labels)
        axs[0,0].set_ylabel(yLabel)
        axs[0,0].yaxis.label.set(fontsize=15)
        axs[0,0].set_xticks(x_line_positions,minor=True)
        axs[0,0].xaxis.grid(True,which='minor')
        axs[0,0].yaxis.grid(True,which='major')
        #------------------------------------------------------------------
#        axs[0,0].set_xticks(tick_positions)
#        axs[0,0].set_xticklabels(tick_labels_double_X, fontsize=label_fontSize)
        #------------------------------------------------------------------
    elif ii == 1:
    #elif ii == 2:
        mesh1 = axs[0,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,1].title.set_text(season_titles[ii])
        axs[0,1].set_yticks(y_ticks)
        axs[0,1].set_yticklabels(y_tick_labels)
        axs[0,1].set_xticks(x_line_positions,minor=True)
        axs[0,1].xaxis.grid(True,which='minor')
        axs[0,1].yaxis.grid(True,which='major')
        #axs[0,1].set_ylabel(yLabel)
        #axs[0,1].yaxis.label.set(fontsize=15)
        #------------------------------------------------------------------
        #------------------------------------------------------------------
    elif ii == 2:
    #elif ii == 3:
        mesh1 = axs[1,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,0].title.set_text(season_titles[ii])
        axs[1,0].set_yticks(y_ticks)
        axs[1,0].set_yticklabels(y_tick_labels)
        axs[1,0].set_ylabel(yLabel)
        axs[1,0].yaxis.label.set(fontsize=15)
        axs[1,0].set_xticks(x_line_positions,minor=True)
        axs[1,0].xaxis.grid(True,which='minor')
        axs[1,0].yaxis.grid(True,which='major')
        #------------------------------------------------------------------
        #------------------------------------------------------------------
    else:
        mesh1 = axs[1,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,1].title.set_text(season_titles[ii])
        axs[1,1].set_yticks(y_ticks)
        axs[1,1].set_yticklabels(y_tick_labels)
        axs[1,1].set_xticks(x_line_positions,minor=True)
        axs[1,1].xaxis.grid(True,which='minor')
        axs[1,1].yaxis.grid(True,which='major')
        #axs[1,1].set_ylabel(yLabel)
        #axs[1,1].yaxis.label.set(fontsize=15)
        #------------------------------------------------------------------
        #------------------------------------------------------------------

#stats_10 = list(range(int(np.shape(stats_vectors)[0])))
#stats_90 = list(range(int(np.shape(stats_vectors)[0])))

cbar_label = "Log base 10 of probability"
cbar_fontSize = 15
cbar_nBins = 20

#cbaxes = fig.add_axes([0.95,0.1, 0.03, 0.8])

#cbar = plt.colorbar(mesh1, ax=axs.ravel().tolist(), cax = cbaxes)
cbar = plt.colorbar(mesh1, ax=axs.ravel().tolist())
cbar.ax.set_ylabel(cbar_label, fontsize = cbar_fontSize)
cbar_tick_labels_pre = [float(t.get_text().replace('−','-')) for t in cbar.ax.get_yticklabels()]
cbar_round = 4
cbar_tick_labels = [round(10**t, cbar_round) for t in cbar_tick_labels_pre]

cbar_label_2 = "probability"

# define functions that relate the two colorbar scales
# e.g., Celcius to Fahrenheit and vice versa
def logP_to_P(x):
    val = 10**(np.ma.masked_invalid(x))
    return np.ma.masked_invalid(val)
def P_to_logP(x):
    val = np.log10(np.ma.masked_invalid(x))
    return np.ma.masked_invalid(val)

# create a second axes
cbar2 = cbar.ax.secondary_yaxis('left',functions=(logP_to_P,P_to_logP))
cbar2.set_ylabel(cbar_label_2, fontsize = cbar_fontSize, y=0.87)
#cbar2.locator_params(nbins=cbar_nBins)
cbar2.yaxis.set_major_locator(LogLocator(subs=(1,.75,.5,.25)))
cbar2.yaxis.set_major_formatter('{x:g}')

fig.suptitle(fig_fullTitle)

plt.show()



