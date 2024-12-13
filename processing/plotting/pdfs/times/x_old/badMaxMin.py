# Copied from "plot_pdf_seasonal_island_swap_v1.py"

#------------------------------------------------------------------
pdf_file_name = "pdf_data_output_seasonal_oneFileTest.p"
#------------------------------------------------------------------



import pickle
import numpy as np
import matplotlib.pyplot as plt
import datetime

#------------------------------------------------------------------
base_path = '/home/blaughli/tracking_project/'
pdf_raw_directory = base_path + 'practice/bounding_boxes/final_locations/z_output/'

pdf_raw_file = pdf_raw_directory + pdf_file_name
#------------------------------------------------------------------

#base_datetime = datetime.datetime(1970,1,1,0,0,0)

file = open(pdf_raw_file,'rb')
pdf_list_connectivity,pdf_list_settleTime,settlement_boxes_test_array,settlement_times_test_array,counter_array = pickle.load(file)  # When the new calc is done, saved a counter_array for checking consistency
file.close()





#pdf_conn_list = []
#for pdf in pdf_list_connectivity:
#    pdf_conn_list.append(np.log10(np.transpose(pdf)))


pdf_max_val = 0
pdf_min_val = 999999
minFindSwitch = 1
#for pdf in pdf_list_settleTime:
# Can just use the full pdf (not seasonal) since it has all of the same info
pdf = pdf_list_settleTime[0]    
for ii in range((np.shape(pdf)[1])):
    if sum(pdf[:,ii]) > 0:
        pdf_max_val = ii
        if minFindSwitch == 1:
            pdf_min_val = ii
            minFindSwitch = 0


print(pdf_max_val)
print(pdf_min_val)

# Determined elsewhere (see/run "check_box_numbers.py")
first_continent_box_dex = 20
num_dummy_lines = 1

tick_positions = [1,5,9,10,13,16,17,19,23,29,32,37,41,47,56,60,67,78]
tick_labels = ['SC','C','SB','SN','SR','A','SC','SM','TJ','PV','PM','PC','PB','CB','PR','PA','CM','CB']

n_boxes_seeded = int(np.shape(pdf_list_settleTime[0])[0])
n_times = int(np.shape(pdf_list_settleTime[0])[1])
X = np.arange(-0.5, n_boxes_seeded, 1)
Y = np.arange(-0.5, n_times, 1)

fig,axs = plt.subplots(2,2)
#plt.setp(axs,xticks=tick_positions,xticklabels=tick_labels,yticks=tick_positions,yticklabels=tick_labels)


for ii in range(1,len(pdf_list_settleTime)):
#for pdf_plot in pdf_list:

    print(ii)

    if ii == 1:
        mesh1 = axs[0,0].pcolormesh(X,Y,pdf_list_settleTime[ii].T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,0].title.set_text("Winter (DJF)")
    elif ii == 2:
        mesh1 = axs[0,1].pcolormesh(X,Y,pdf_list_settleTime[ii].T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,1].title.set_text("Spring (MAM)")
    elif ii == 3:
        mesh1 = axs[1,0].pcolormesh(X,Y,pdf_list_settleTime[ii].T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,0].title.set_text("Summer (JJA)")
    else:
        mesh1 = axs[1,1].pcolormesh(X,Y,pdf_list_settleTime[ii].T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,1].title.set_text("Fall (SON)")

    plt.colorbar(mesh1)

fig.suptitle("PDFs of Settlement Time (y-axis) vs Release Location (x-axis).\n Grouped according to season of release")

plt.show()



