# Copied from "stablePlotter_update2", modified to handle multiple input files from the parallelized connectivity prePdf calculation

script_version = "allPLDs_v1"



# ----------------------------------------------
# ----------------------------------------------
# Only for testing
# ----------------------------------------------
# ----------------------------------------------

# ie these should be input arguments

#logSwitch = False
logSwitch = True

annualOnlySwitch = False
#annualOnlySwitch = True

#pld_index = 5
#pld_index = 0

ignoreStagnantSwitch = False
#ignoreStagnantSwitch = True

plotDiagonalLineSwitch = False
#plotDiagonalLineSwitch = True

annualSeasonIndex = 4

tracking_dir = "/data03/blaughli/tracking_output/Mercator_coastalCells_1993_2018_kickSTD_0p0___global-reanalysis-phy-001-030-daily_1993_2018"
#tracking_dir = "/data03/blaughli/tracking_output/z_complete/z_tests/Mercator_coastalCells_1993_1993_kick_0p0_MaskDepthInForcing___Mercator"
#tracking_dir = "/home/blaughli/tracking_project_v2/t_scraps/test_dir_con"


subplot_label_fontsize = 6
title_fontsize = 8


#connectivity_data_dir = args.connectivitydatadir
#pDrake_switch = 

# New approach: set min val to log(0.0001), so 0.0001 is our smallest colorbar tick. indicate that in colorbar tick labels
pdf_min_val = 0
pdf_min_val_log = 0.0001

pdf_max_val = 0.01
pdf_max_val_log = 0.1


# ----------------------------------------------
# ----------------------------------------------
# ----------------------------------------------
# ----------------------------------------------


# Also set vmax to log(0.1), to match Patricks's 2011 figures
#pdf_max_val = 0.02
#pdf_max_val = 0.01
#pdf_max_val = 0.1
#pdf_max_val = 0.2
#---------------------------------------------------------------------

import matplotlib.colors as colors
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import os
from pathlib import Path
from glob import glob

#parser = argparse.ArgumentParser()
#parser.add_argument("trackingdir")
###parser.add_argument("connectivitydatadir")
#parser.add_argument("pdrakeswitch", nargs='?', type=str)
#args = parser.parse_args()

#tracking_dir = args.trackingdir
###connectivity_data_dir = args.connectivitydatadir
#pDrake_switch = args.pdrakeswitch

input_connectivity_dir_stem = "v_connHist_files"
output_csv_dir_stem = "y_csv_files"
output_figures_dir_stem = "z_figures"

connectivity_data_dir = os.path.join(tracking_dir,input_connectivity_dir_stem)


output_figures_dir = os.path.join(connectivity_data_dir,output_figures_dir_stem)
Path(output_figures_dir).mkdir(parents=True, exist_ok=True)

output_csv_dir = os.path.join(connectivity_data_dir,output_csv_dir_stem)
Path(output_csv_dir).mkdir(parents=True, exist_ok=True)

figure_title_base = Path(tracking_dir).stem

connectivity_pdf_list = []
release_counts_per_polygon_array_list = []

pdf_max_val = -999999
pdf_min_val = 999999
release_months_list_pre = np.array(["Dec-Feb", "March-May","June-Aug","Sep-Nov","Jan-Dec"])
#release_months_list_pre = ["Dec-Feb", "March-May","June-Aug","Sep-Nov","Jan-Dec"]

if annualOnlySwitch:
    pdfs_to_plot_indices = [4]
else:
    pdfs_to_plot_indices = [0,1,2,3]

release_months_list = release_months_list_pre[pdfs_to_plot_indices]
    

# HACK FOR TESTING ONLY
#pld_to_plot = 0

num_seasons = len(pdfs_to_plot_indices)


num_released_per_season = []
num_settled_per_season = []
for season_index in range(num_seasons):
    num_released_per_season.append(0)
    num_settled_per_season.append(0)

num_released_per_season = np.array(num_released_per_season)
num_settled_per_season = np.array(num_settled_per_season)

file_list = sorted(glob(os.path.join(connectivity_data_dir,"*.npz")))

d = np.load(file_list[0])
pld_list = d['pld_list']

for pld_index in range(len(pld_list)):
    #pld_string = f"{pld_list[pld_index][0]}_{pld_list[pld_index][1]}"
    pld_string = f"{pld_list[pld_index][0]:03}_{pld_list[pld_index][1]:03}"
    
    print(f"PLD: {pld_string}")

    file_counter = 0

    #for connectivity_data_file_byte in os.listdir(directory):
    for connectivity_data_file in file_list:

#        if file_counter > 1:
#            break


    #        print('hi4')

        file_counter += 1


    #        print(f"FILE: {connectivity_data_file}")
    #        print(f"file counter: {file_counter}")
        print(f"file {file_counter:>3}/{len(file_list)}")
            
        d = np.load(connectivity_data_file)
        ###pld = d['pld_days']
        release_counts_per_polygon_array = d['release_counts_per_polygon_array']
        if ignoreStagnantSwitch:
            prePdf_arrays_connectivity = d['prePdf_arrays_connectivity_noStagnation']
        else:
            prePdf_arrays_connectivity = d['prePdf_arrays_connectivity']

        if file_counter == 1:
            
            # Setup stuff, now that we're combining multiple prePdfs 
        #    d = np.load(connectivity_data_file)
            #pld_list = d['pld_list']
            ###pld = d['pld_days']
        #    release_counts_per_polygon_array = d['release_counts_per_polygon_array']

        #    if ignoreStagnantSwitch:
        #        prePdf_arrays_connectivity = d['prePdf_arrays_connectivity_noStagnation']
        #    else:
        #        prePdf_arrays_connectivity = d['prePdf_arrays_connectivity']

            


            connectivity_file_stem = Path(connectivity_data_file).stem

            # The figure titles also need to change when I start plotting seasons again.  So, this is bugged
            # for seasons

            if annualOnlySwitch:
                figure_file_leaf = f"{figure_title_base}_pld_{pld_string}_annual.png"
                #figure_file_leaf = figure_title_base + "_annual" + '.png'
            else:
                figure_file_leaf = f"{figure_title_base}_pld_{pld_string}_seasonal.png"
                #figure_file_leaf = figure_title_base + "_seasons" + '.png'
            figure_file = os.path.join(output_figures_dir,figure_file_leaf)

            csv_file_leaf = connectivity_file_stem + '.png'
            csv_file = os.path.join(output_csv_dir,csv_file_leaf)

    #            wordcount_line_limit = 75
    #            num_lines_filename = int(np.ceil(len(connectivity_file_stem)/wordcount_line_limit))
    #            figure_title = ""
            
            ### This code will not really work for seasonal plots... been debugging with the annual pdf for so long, I've
            ### written bugs into this script.

    #            for ii in range(num_lines_filename):
    #                figure_title += f"{figure_title_base[ii*wordcount_line_limit:(ii+1)*wordcount_line_limit]}"
    #                if ii < num_lines_filename - 1:
    #                    figure_title += " ...\n"

        #release_counts_per_polygon_array_list = []
    #           print('hi2')
        
            if annualOnlySwitch:
                pdf = np.copy(prePdf_arrays_connectivity[annualSeasonIndex,pld_index])
                connectivity_pdf_list.append(np.zeros_like(pdf))
                release_counts_per_polygon_array_list.append(np.zeros_like(release_counts_per_polygon_array[annualSeasonIndex]))
            else: 
                for season_index in pdfs_to_plot_indices:
                    pdf = np.copy(prePdf_arrays_connectivity[season_index,pld_index])
                    connectivity_pdf_list.append(np.zeros_like(pdf))
                    release_counts_per_polygon_array_list.append(np.zeros_like(release_counts_per_polygon_array[season_index]))
          

        if annualOnlySwitch:
            pdf = np.copy(prePdf_arrays_connectivity[annualSeasonIndex,pld_index])
            connectivity_pdf_list[0] += pdf
            release_counts_per_polygon_array_list[0] += release_counts_per_polygon_array[annualSeasonIndex]
            
            num_released_per_season[0] += int(np.sum(release_counts_per_polygon_array[annualSeasonIndex])) 
            num_settled_per_season[0] += int(np.sum(prePdf_arrays_connectivity[annualSeasonIndex,pld_index]))
        
        else: 
            for season_index in pdfs_to_plot_indices:
                pdf = np.copy(prePdf_arrays_connectivity[season_index,pld_index])
                connectivity_pdf_list[season_index] += pdf
                release_counts_per_polygon_array_list[season_index] += release_counts_per_polygon_array[season_index]
            
                num_released_per_season[season_index] += int(np.sum(release_counts_per_polygon_array[season_index])) 
                num_settled_per_season[season_index] += int(np.sum(prePdf_arrays_connectivity[season_index,pld_index]))
            



    for adjusted_season_index in pdfs_to_plot_indices:
        # clunky shifting to account for the annualOnly case
        if annualOnlySwitch:
            adjusted_season_index = 0

    #    adjusted_season_index = 0 #DJF
    #    adjusted_season_index = 1 #MAM
    #    adjusted_season_index = 2 #JJA
    #    adjusted_season_index = 3 #SON

        # Why did I need to access element 0 of <release_counts_per_polygon_array_list[adjusted_season_index]> before?  It seems like I've gone back and forth on this, need to understand

    #    if annualOnlySwitch:
    #        connectivity_pdf_list[adjusted_season_index] = connectivity_pdf_list[adjusted_season_index] / release_counts_per_polygon_array_list[adjusted_season_index][:, np.newaxis] 
    #        ####connectivity_pdf_list[0] = connectivity_pdf_list[0] / release_counts_per_polygon_array_list[0][0][:, np.newaxis] # Now need to access element 0
    #        # Testing
    #        #connectivity_pdf_list[adjusted_season_index] = connectivity_pdf_list[adjusted_season_index] 
    #    else:
    #        connectivity_pdf_list[adjusted_season_index] = connectivity_pdf_list[adjusted_season_index] / release_counts_per_polygon_array_list[adjusted_season_index][:, np.newaxis] 
    #        #connectivity_pdf_list[adjusted_season_index] = connectivity_pdf_list[adjusted_season_index] / release_counts_per_polygon_array_list[adjusted_season_index][0][:, np.newaxis] # Now need to access element 0
    #        # Testing
    #        #connectivity_pdf_list[adjusted_season_index] = connectivity_pdf_list[adjusted_season_index] 
            
        connectivity_pdf_list[adjusted_season_index] = connectivity_pdf_list[adjusted_season_index] / release_counts_per_polygon_array_list[adjusted_season_index][:, np.newaxis] 

        if np.amax(connectivity_pdf_list[adjusted_season_index]) > pdf_max_val:
            pdf_max_val = np.amax(connectivity_pdf_list[adjusted_season_index])
        if np.amin(np.ma.masked_invalid(connectivity_pdf_list[adjusted_season_index])) < pdf_min_val:
            pdf_min_val = np.amin(np.ma.masked_invalid(connectivity_pdf_list[adjusted_season_index]))




    settle_strength_array = num_settled_per_season/num_released_per_season

    # This should be the same as the settle strength for the AnnualOnly option!!!  CHECK!!!
    settle_strength_full_run = np.sum(num_settled_per_season)/np.sum(num_released_per_season)



    if annualOnlySwitch:
        nrows = 1
        ncols = 1
        #fig,axs = plt.subplots(1,1, squeeze=False, figsize=fig_size)
        #fig,axs = plt.subplots(1,1)
        #fig,axs = plt.subplots(1,1, figsize = fig_size)
    else:
        #fig,axs = plt.subplots(2,2, squeeze=False, sharex=True, sharey=True, figsize=fig_size)
        nrows = 2
        ncols = 2
        #fig,axs = plt.subplots(2,2, squeeze=False)
        #fig,axs = plt.subplots(2,2)
        #fig,axs = plt.subplots(2,2, figsize = fig_size)


    label_fontsize=6
    #fig_size = (16,9)
    #fig_size = (6*ncols+2,6*nrows)
    fig_size = (6*ncols+2,6*nrows)


    ###fig,axs = plt.subplots(nrows=nrows,ncols=ncols, squeeze=False, sharex=True, sharey=True, figsize=fig_size)
    #fig,axs = plt.subplots(nrows=nrows,ncols=ncols, squeeze=False, sharex=True, sharey=True, figsize=fig_size, gridspec_kw={'hspace':0.1, 'wspace':0.1})
    #fig,axs = plt.subplots(nrows=nrows,ncols=ncols, squeeze=False, figsize=fig_size, gridspec_kw={'hspace':0.1, 'wspace':0.1})
    fig,axs = plt.subplots(nrows=nrows,ncols=ncols, squeeze=False, figsize=fig_size, gridspec_kw={'hspace':0.0001, 'wspace':0.15})

    pcolormesh_plot_list = []


    num_dummy_lines = 1

    # Wait, there are 5 pdfs in the list - the first (index 0) is the overall pdf (non-seasonal).  So, do I have a 1-off error here?




    boundary_index = 27
    island_index = 489

    print(f"Number of pdfs: {len(connectivity_pdf_list)}")

    for season_index in pdfs_to_plot_indices:
        
        # What was I doing here?  
    #    if season_index == num_seasons -1:
    #        axis_index = 0
    #    else:
    #        axis_index = season_index

        if annualOnlySwitch:
            season_index = 0


        pdf_plot = connectivity_pdf_list[season_index]
        #pdf_plot = connectivity_pdf_list[season_index][0,:,:]   # Why is this neccesary now???

        #subplot_title = ""
        #subplot_title += f"PLD: {pld_list[pld_index][0]}-{pld_list[pld_index][1]} days"
        #subplot_title += f"\nReleased {release_months_list[season_index]}, 0-20m"
        #subplot_title += f"\nnum floats: {num_released_per_season[season_index]}, n settlers: {num_settled_per_season[season_index]}, settle strength: {settle_strength_array[season_index]:.3f}"
        
        subplot_label = f"num floats: {num_released_per_season[season_index]}, n settlers: {num_settled_per_season[season_index]}, settle strength: {settle_strength_array[season_index]:.3f}"
        #subplot_label = f"\nnflts: {num_released_per_season[season_index]}, n settlers: {num_settled_per_season[season_index]}, settle strength: {settle_strength_array[season_index]:.3f}"

    #if np.logical_not(annualOnlySwitch):
    #        subplot_title += f"\nn Overall: floats: {np.sum(num_released_per_season)}, n settlers: {np.sum(num_settled_per_season)}, settle strength: {settle_strength_full_run:.3f}"


        n_boxes_seeded = int(np.shape(pdf_plot)[1])
        n_boxes_settled = int(np.shape(pdf_plot)[0])
        X = np.arange(-0.5, n_boxes_settled, 1)
        Y = np.arange(-0.5, n_boxes_seeded, 1)

        #if pDrake_switch is None:
        #    pdf_plot[0:boundary_index,:] = 0
        #    pdf_plot[:,0:boundary_index] = 0
        #    pdf_plot[island_index - boundary_index:island_index,:] = 0
        #    pdf_plot[:,island_index - boundary_index:island_index] = 0

        if season_index == 0:
            conn_csv_file_name = output_csv_dir + "/" + connectivity_file_stem + f"_annual.csv"
            #conn_csv_file_name = output_csv_dir_stem + "/" + connectivity_file_stem + f"_annual.csv"
        else:
            conn_csv_file_name = output_csv_dir + "/" + connectivity_file_stem + f"_season_{season_index + 1}/4.csv"
            #conn_csv_file_name = output_csv_dir_stem + "/" + connectivity_file_stem + f"_season_{season_index + 1}/4.csv"

        ###np.savetxt('connectivity.csv', pdf_plot, delimiter = ',')
        #np.savetxt(conn_csv_file_name, pdf_plot, delimiter = ',')

        ###pcolormesh_plot_list.append(ax.pcolormesh(X,Y,pdf_plot.T,cmap='jet',norm=colors.LogNorm(vmin=0.001,vmax=.1),shading='auto'))
        #pcolormesh_plot_list.append(ax[axis_index].pcolormesh(X,Y,pdf_plot.T,cmap='jet',norm=colors.LogNorm(vmin=0.001,vmax=.1),shading='auto'))
        
        ###ax.plot([axis_index,np.shape(pdf_plot)[1]-1],[0,np.shape(pdf_plot)[0]-1],color="white",linewidth=0.5)
        #ax.plot([axis_index,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)

        if annualOnlySwitch:
            if logSwitch:
                mesh1 = axs[0,0].pcolormesh(X,Y,pdf_plot.T,cmap='jet',norm=colors.LogNorm(vmin=pdf_min_val_log,vmax=pdf_max_val_log),shading='auto')
            else:
                mesh1 = axs[0,0].pcolormesh(X,Y,pdf_plot.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val,shading='auto')
            if plotDiagonalLineSwitch:
                #axs[0,0].plot([axis_index,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)
                axs[0,0].plot([season_index,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)
            axs[0,0].set_title("Annual data")
            #axs[0,0].set(xlabel = "Annual data")
            #axs[0,0].set_aspect('equal')
            axs[0,0].set_aspect(1.0)
            break

        if season_index == 0:
            if logSwitch:
                mesh1 = axs[0,0].pcolormesh(X,Y,pdf_plot.T,cmap='jet',norm=colors.LogNorm(vmin=pdf_min_val_log,vmax=pdf_max_val_log),shading='auto')
            else:
                mesh1 = axs[0,0].pcolormesh(X,Y,pdf_plot.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val,shading='auto')
            if plotDiagonalLineSwitch:
                #axs[0,0].plot([axis_index,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)
                axs[0,0].plot([season_index,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)
            #axs[0,0].title.set_text("Winter (DJF)")
            #axs[0,0].set_xlabel(f"Winter (DJF)  {subplot_label}", size=subplot_label_fontsize)
            axs[0,0].set_title(f"Winter (DJF)  {subplot_label}", size=subplot_label_fontsize)
            #axs[0,0].set_aspect('equal')
            axs[0,0].set_aspect(1.0)
        elif season_index == 1:
            if logSwitch:
                mesh2 = axs[0,1].pcolormesh(X,Y,pdf_plot.T,cmap='jet',norm=colors.LogNorm(vmin=pdf_min_val_log,vmax=pdf_max_val_log),shading='auto')
            else:
                mesh2 = axs[0,1].pcolormesh(X,Y,pdf_plot.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val,shading='auto')
            if plotDiagonalLineSwitch:
                axs[0,1].plot([season_index,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)
                #axs[0,1].plot([axis_index,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)
            #axs[0,1].title.set_text("Spring (MAM)")
            #axs[0,1].set_xlabel(f"Spring (MAM)  {subplot_label}", size=subplot_label_fontsize)
            axs[0,1].set_title(f"Spring (MAM)  {subplot_label}", size=subplot_label_fontsize)
            #axs[0,1].set_aspect('equal')
            axs[0,1].set_aspect(1.0)
        elif season_index == 2:
            if logSwitch:
                mesh3 = axs[1,0].pcolormesh(X,Y,pdf_plot.T,cmap='jet',norm=colors.LogNorm(vmin=pdf_min_val_log,vmax=pdf_max_val_log),shading='auto')
            else:
                mesh3 = axs[1,0].pcolormesh(X,Y,pdf_plot.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val,shading='auto')
            if plotDiagonalLineSwitch:
                axs[1,0].plot([season_index,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)
                #axs[1,0].plot([axis_index,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)
            #axs[1,0].title.set_text("Summer (JJA)")
            #axs[1,0].set_xlabel(f"Summer (JJA)  {subplot_label}", size=subplot_label_fontsize)
            axs[1,0].set_title(f"Summer (JJA)  {subplot_label}", size=subplot_label_fontsize)
            #axs[1,0].set_aspect('equal')
            axs[1,0].set_aspect(1.0)
        else:
            if logSwitch:
                mesh4 = axs[1,1].pcolormesh(X,Y,pdf_plot.T,cmap='jet',norm=colors.LogNorm(vmin=pdf_min_val_log,vmax=pdf_max_val_log),shading='auto')
            else:
                mesh4 = axs[1,1].pcolormesh(X,Y,pdf_plot.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val,shading='auto')
            if plotDiagonalLineSwitch:
                axs[1,1].plot([axis_index,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)
            #axs[1,1].title.set_text("Fall (SON)")
            #axs[1,1].set_xlabel(f"Fall (SON)  {subplot_label}", size=subplot_label_fontsize)
            axs[1,1].set_title(f"Fall (SON)  {subplot_label}", size=subplot_label_fontsize)
            #axs[1,1].set_aspect('equal')
            axs[1,1].set_aspect(1.0)






        #axs.set_title(subplot_title)
        #axs.title.set(subplot_title)
        #ax.title.set(subplot_title)

        #if annualOnlySwitch:
        #    break


    overall_plot_title = ""
    overall_plot_title += f"PLD: {pld_list[pld_index][0]}-{pld_list[pld_index][1]} days"
    overall_plot_title += f"\nReleased {release_months_list_pre[-1]}, 0-20m"
    #overall_plot_title += f"\nReleased {release_months_list[season_index]}, 0-20m"
    #overall_plot_title += f"\nnum floats: {num_released_per_season[season_index]}, n settlers: {num_settled_per_season[season_index]}, settle strength: {settle_strength_array[season_index]:.3f}"
    overall_plot_title += f"\nnum floats: {np.sum(num_released_per_season)}, n settlers: {np.sum(num_settled_per_season)}, settle strength: {settle_strength_full_run:.3f}"


    cbar_label = "probability"

    cbar = plt.colorbar(mesh1, ax=axs.ravel(), shrink=0.6)
    #cbar = plt.colorbar(mesh1, ax=axs.ravel(), shrink=0.8)
    #cbar = plt.colorbar(mesh1, ax=axs.ravel())
    #cbar = plt.colorbar(mesh1, ax=axs.ravel().tolist())

    #cbar = plt.colorbar(pcolormesh_plot_list[0],label=cbar_label,extend="both")

    #axs.set_aspect('equal')

    figure_title = figure_title_base + "\n" + overall_plot_title

    #fig.suptitle(figure_title, fontsize=title_fontsize)
    fig.suptitle(figure_title, fontsize=title_fontsize, y=0.92)

    #plt.figtext(0.5, 0.02, figure_title_base, ha="center")
    ###plt.figtext(0.5, 0.02, figure_title, ha="center")

    plt.savefig(figure_file, bbox_inches = "tight")


    #plt.show()



