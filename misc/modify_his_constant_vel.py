
import netCDF4
import os
import numpy as np
import time


#file_directory = '/home/blaughli/tracking_project_v2/misc/dummy_files/zero_vel_10days'
file_directory = '/home/blaughli/forcing_samples/onshore_constant_vel_10days_1990'
#file_directory = '/home/blaughli/forcing_samples/single_year_1995/'

constant_velocity_value = 1

# Set the number of timesteps of the netCDF file to read at once, since we can't read it all at once bc of memory limitations
n_times_read = 10

vars_to_modify = ['u','ubar','v','vbar']
#vars_to_modify = ['u','ubar','v','vbar','w']


for filename_pre in os.listdir(file_directory):
    
    if filename_pre.endswith('.nc'):

        filename = os.path.join(file_directory,filename_pre)

        print(filename)

        # Maybe if I do this one variable at a time, I won't blow out the memory..?

        dset = netCDF4.Dataset(filename, 'r+')

        t_start_var = time.time()

        var_dex = 0

        for var in vars_to_modify:

            var_dex += 1

            time_dimension_length = np.shape(dset[var])[0]
            
            for chunk_start in range(0,time_dimension_length,n_times_read):
    
                t_start_chunk = time.time()

                #var_orig = dset[var][:] 
                #var_mask = var_orig.mask
                #var_data = var_orig.data
                
                #var_data_new = np.where(var_mask == False, constant_velocity_value, var_data)
                #var_new = np.ma.masked_array(var_data_new, var_mask)
               
                if chunk_start + n_times_read >= time_dimension_length:
                    #print(f"{chunk_start}, {chunk_start + n_times_read}, {time_dimension_length}    last step")
                    if var[-3:] == 'bar':
                        dset[var][chunk_start:,:,:] = np.ma.masked_array(np.where(dset[var][chunk_start:,:,:].mask == False, constant_velocity_value, dset[var][chunk_start:,:,:].data), dset[var][chunk_start:,:,:].mask)
                    else:
                        dset[var][chunk_start:,:,:,:] = np.ma.masked_array(np.where(dset[var][chunk_start:,:,:,:].mask == False, constant_velocity_value, dset[var][chunk_start:,:,:,:].data), dset[var][chunk_start:,:,:,:].mask)
                else:
                    #print(f"{chunk_start}, {chunk_start + n_times_read}, {time_dimension_length}")
                    if var[-3:] == 'bar':
                        dset[var][chunk_start:chunk_start+n_times_read,:,:] = np.ma.masked_array(np.where(dset[var][chunk_start:chunk_start+n_times_read,:,:].mask == False, constant_velocity_value, dset[var][chunk_start:chunk_start+n_times_read,:,:].data), dset[var][chunk_start:chunk_start+n_times_read,:,:].mask)
                    else:
                        dset[var][chunk_start:chunk_start+n_times_read,:,:,:] = np.ma.masked_array(np.where(dset[var][chunk_start:chunk_start+n_times_read,:,:,:].mask == False, constant_velocity_value, dset[var][chunk_start:chunk_start+n_times_read,:,:,:].data), dset[var][chunk_start:chunk_start+n_times_read,:,:,:].mask)
                


                t_end = time.time()
                chunk_runtime = t_end - t_star_chunkt
                chunk_runtime_minutes = chunk_runtime/60
                #chunk_runtime_minutes = int(np.round(report_runtime/60))

                print(f"current var chunk time: {chunk_runtime_minutes:02.2f} minutes")
           
            print(f"{var_dex}/{len(vars_to_modify)} variables done")
            chunked_var_runtime = t_end - t_start_var
            chunked_var_runtime_minutes = chunked_var_runtime/60
            print(f"current chunked var time: {chunked_var_runtime_minutes:02.2f} minutes")




                #dset[var][:] = np.ma.masked_array(np.where(dset[var][:].mask == False, constant_velocity_value, dset[var][:].data), dset[var][:].mask)
                #dset[var][:] = np.ma.masked_array(np.where(var_orig.mask == False, constant_velocity_value, var_orig.data), var_orig.mask)
               
            #break

    #    u_orig = dset['u'][:] 
    #    u_mask = u_orig.mask
    #    u_data = u_orig.data
    #    u_data_new = np.where(u_mask == False, constant_velocity_value, u_data)
    #    u_new = np.ma.masked_array(u_data_new, u_mask)
    #    
    #    ubar_orig = dset['ubar'][:] 
    #    ubar_mask = ubar_orig.mask
    #    ubar_data = ubar_orig.data
    #    ubar_data_new = np.where(ubar_mask == False, constant_velocity_value, ubar_data)
    #    ubar_new = np.ma.masked_array(ubar_data_new, ubar_mask)
    #
    #    v_orig = dset['v'][:] 
    #    v_mask = v_orig.mask
    #    v_data = v_orig.data
    #    v_data_new = np.where(v_mask == False, constant_velocity_value, v_data)
    #    v_new = np.ma.masked_array(v_data_new, v_mask)
    #
    #    vbar_orig = dset['vbar'][:] 
    #    vbar_mask = vbar_orig.mask
    #    vbar_data = vbar_orig.data
    #    vbar_data_new = np.where(vbar_mask == False, constant_velocity_value, vbar_data)
    #    vbar_new = np.ma.masked_array(vbar_data_new, vbar_mask)
    #
        #dset['w'][:] = constant_velocity_value

        dset.close()










