% to call from command line:
%   matlab -nodisplay -nosplash -nodesktop -r "run('s02_ONLY_RUN_ONCE_calc_dist2coast_ISLANDS_Mercator.m');exit;"

% SYNTAX
% ------
% dist2cst=dist2coast_gcircle(lonr,latr,maskr,lonc,latc);
%
% where lonr,latr are the lon/lat points of the grid 
%       maskr is the mask of the grid (used to only find dist for ocean points
%       lonc,latc define the coastline (can and should probably have NaNs in it)
%  and  dist is what's returned (in meters);

% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
% --------------------------- EDIT THESE --------------------------------------------
% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
addpath('/home/blaughli/tracking_project_v2/misc/z_boxes/y_external_matlab_scripts')


grid_file = '/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid_Matlab.mat'
load(grid_file)

islands_input_dir = '/home/blaughli/tracking_project_v2/misc/z_boxes/a_islands/z_output/'

output_dir = '/home/blaughli/tracking_project_v2/misc/z_boxes/a_islands/z_output/'

num_islands = 3

for island_dex = 1:num_islands

    save_file = strcat('dist_2_coast_field_coastline_Mercator_island_number_',num2str(island_dex),'.mat')
    save_file_path = strcat(output_dir,save_file)

    coast_lon_file = strcat('coastline_coords_Mercator_lon_island_number_',num2str(island_dex),'.txt')
    coast_lat_file = strcat('coastline_coords_Mercator_lat_island_number_',num2str(island_dex),'.txt')
    coast_lon_in = strcat(islands_input_dir,coast_lon_file)
    coast_lat_in = strcat(islands_input_dir,coast_lat_file)


    % Dirty Matlab - need to add 1 to each element of these coast i/j vectors
    fileID = fopen(coast_lon_in,'r');
    formatSpec = '%f';
    coast_lon = fscanf(fileID,formatSpec);
    fclose(fileID);

    fileID = fopen(coast_lat_in,'r');
    formatSpec = '%f';
    coast_lat = fscanf(fileID,formatSpec);
    fclose(fileID);

    dist_field = dist2coast_gcircle(lon_rho,lat_rho,mask_rho,coast_lon,coast_lat);

    save(save_file_path,'dist_field')

end



