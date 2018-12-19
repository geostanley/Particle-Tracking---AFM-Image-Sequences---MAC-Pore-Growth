%%%=== MACanalysis_TrackVideos_concatonate ===%%%

% This is script 4/5.

% This script concatonates the data from several imaging experiments, and
% saves them all into one data structure. This can then be loaded up into
% the final script for calculating the time of individual proe growths, and
% plotting histograms etc.

%% Input data directory, file name, and output directory

clear variables
close all
clc

%%% Enter all data structure directories and filenames (as strings in cells)
%%% to be concatonated.
Load_directories = [{'Z:\Users\George\Documents\PhD\MAC_Manuscript\MAC_Analysis_Everything\EPFL1_experiment\DataStructures'},...
    {'Z:\Users\George\Documents\PhD\MAC_Manuscript\MAC_Analysis_Everything\EPFL2_experiment\DataStructures'}];

DataStructure_LoadNames = [{'MAC_EPFL1_pore_track_vids_TrueTracks'},...
    {'MAC_EPFL2_pore_track_vids_TrueTracks'}];

% Save concatonated data structure out to:
outdirectory   = 'Z:\Users\George\Documents\PhD\MAC_Manuscript\MAC_Analysis_Everything\Concatonated_data_sets\DataStructures';
DataStructure_SaveName  = 'MAC_EPFL1_EPFL2_concatonated';

%% Create cell of full file names for loading

loadfilename_cell = cell(size(DataStructure_LoadNames));

for i = 1:length(DataStructure_LoadNames)
    
    directory = Load_directories{i};
    name      = DataStructure_LoadNames{i};
    filename  = strcat(name, '.mat');
    
    loadfilename         = fullfile(directory, filename);
    loadfilename_cell{i} = loadfilename;
    
end

%% Concatonate data structures

img_seq_cell_cat        = cell(0,0);
mh_vec_cell_cat         = cell(0,0);
first_track_s_array_cat = zeros(0,0);
template_nm_array_cat   = zeros(0,0);

for i = 1:length(loadfilename_cell)
    
    load(loadfilename_cell{i});
    
    img_seq_cell         = MAC_tracking_data_individual_pore_growth_videos_True_tracks.Cropped_pore_growth_cell_true_tracks;
    mh_vec_cell          = MAC_tracking_data_individual_pore_growth_videos_True_tracks.mean_height_true_tracks_cell;
    first_track_s_array  = MAC_tracking_data_individual_pore_growth_videos_True_tracks.firsttrack_sec_array;
    CropSize_template_nm = MAC_tracking_data_individual_pore_growth_videos_True_tracks.CropSize_template_nm;
    
    % create an array of template sizes for plotting in the final script
    template_array = ones(size(img_seq_cell));
    template_array_nm = template_array .* CropSize_template_nm;
    
    img_seq_cell_cat        = horzcat(img_seq_cell_cat, img_seq_cell);
    mh_vec_cell_cat         = horzcat(mh_vec_cell_cat, mh_vec_cell);
    first_track_s_array_cat = horzcat(first_track_s_array_cat, first_track_s_array);
    template_nm_array_cat   = horzcat(template_nm_array_cat, template_array_nm);
    
    clear img_seq_cell
    clear mh_vec_cell
    clear first_track_s_array
    clear CropSize_template_nm
    
end

%% Save concatonated data structure

display('Saving concatonated data structure, ready for fitting...')

MAC_True_tracks_concatonated.Cropped_pore_growth_cell_true_tracks = img_seq_cell_cat;
MAC_True_tracks_concatonated.mean_height_true_tracks_cell         = mh_vec_cell_cat;
MAC_True_tracks_concatonated.firsttrack_sec_array                 = first_track_s_array_cat;
MAC_True_tracks_concatonated.template_nm_array                    = template_nm_array_cat;

SaveFullFileName = fullfile(outdirectory, strcat(DataStructure_SaveName, '.mat'));
save(SaveFullFileName, 'MAC_True_tracks_concatonated');

% Load this into: MACanalysis_TrackVideos_tau_fitting





