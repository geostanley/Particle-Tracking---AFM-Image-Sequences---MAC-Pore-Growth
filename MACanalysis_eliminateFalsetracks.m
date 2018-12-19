%%%=== MACanalysis_eliminateFalsetracks ===%%%

% This is script 3/5.

% This script loads the individual tracks from one experiment, but only
% saves out the ones designated by the user. This is to eliminate false 
% positive tracks before calculating the time for individual pore formation.

%% Input data directory, file name, and output directory

clear variables
close all
clc

Load_directory = 'Z:\Users\George\Documents\PhD\MAC_Manuscript\MAC_Analysis_Everything\EPFL2_experiment\DataStructures';
outdirectory   = 'Z:\Users\George\Documents\PhD\MAC_Manuscript\MAC_Analysis_Everything\EPFL2_experiment\DataStructures';

DataStructure_LoadName = 'MAC_EPFL2_pore_track_vids';
DataStructure_SaveName = 'MAC_EPFL2_pore_track_vids_TrueTracks';

%%% After reviewing the videos created in the previous script, enter here
%%% the numbers (by index, this is indicated at the end of each file name)
%%% of the pore growths to be kept. I.e., only keep those that you are
%%% certain are the tracks of just one pore, eliminating those for which the
%%% tracking software was incorrect. This includes complete pores that have
%%% drifted into the image due to scanner hysteresis, which the
%%% tracking software cannot distinguish from real pore growing events.

% EPFL 1 data
% True_tracks_idx = [1, 2, 3, 4, 5, 6, 7, 10, 15, 18, 19, 25, 30, 31, 34, 42, 60];

% EPFL 2 data
True_tracks_idx = [1, 4, 6, 7, 8, 9, 10, 11, 13, 14, 20, 21, 23];

%% Load the data structure

display('Loading tracking data structure...')
LoadFullFileName = fullfile(Load_directory, DataStructure_LoadName);
load(strcat(LoadFullFileName, '.mat'));

FrameRate_s              = MAC_tracking_data_individual_pore_growth_videos.FrameRate_s;
CropSize_template_nm     = MAC_tracking_data_individual_pore_growth_videos.CropSize_template_nm;
mean_height_cell         = MAC_tracking_data_individual_pore_growth_videos.mean_height_cell;
Cropped_pore_growth_cell = MAC_tracking_data_individual_pore_growth_videos.Cropped_pore_growth_cell;
trackingtable_cell       = MAC_tracking_data_individual_pore_growth_videos.trackingtable_cell;

%% Create new cell arrays containing information only on true pore growth events

firsttrack_framenumb_array = zeros(size(True_tracks_idx));
firsttrack_sec_array       = zeros(size(True_tracks_idx));

Cropped_pore_growth_cell_true_tracks = cell(size(True_tracks_idx));
mean_height_true_tracks_cell         = cell(size(True_tracks_idx));
    
for i = 1:length(True_tracks_idx)

    track_numb = True_tracks_idx(i);       

    % save out corresponding movie cell
    Cropped_pore_growth_cell_true_tracks{i} = Cropped_pore_growth_cell{track_numb};

    track_info           = trackingtable_cell{track_numb};
    firsttrack_framenumb = track_info(1, 3);
    firsttrack_sec       = firsttrack_framenumb * FrameRate_s;

    firsttrack_framenumb_array(i) = firsttrack_framenumb;
    firsttrack_sec_array(i)       = firsttrack_sec;

    mean_height_array = mean_height_cell{track_numb};
    time_vector       = mean_height_array(:,1);
    ave_height_vector = mean_height_array(:,2);

    % each array is an Nx2 array. Column 1 is time (s) and column 2 is mean
    % frame height (nm)
    mean_height_true_tracks_cell{i} = [time_vector, ave_height_vector];

end

%% Save out the pore growth image sequenecs with mean height arrays and time vectors
% these will be used in the final script to calcualate the time for each
% pore growth event

display('Saving final pore growth video sequences, ready for fitting...')

MAC_tracking_data_individual_pore_growth_videos_True_tracks.Cropped_pore_growth_cell_true_tracks = Cropped_pore_growth_cell_true_tracks;
MAC_tracking_data_individual_pore_growth_videos_True_tracks.mean_height_true_tracks_cell         = mean_height_true_tracks_cell;
MAC_tracking_data_individual_pore_growth_videos_True_tracks.firsttrack_sec_array                 = firsttrack_sec_array;
MAC_tracking_data_individual_pore_growth_videos_True_tracks.CropSize_template_nm                 = CropSize_template_nm;

SaveFullFileName = fullfile(outdirectory, strcat(DataStructure_SaveName, '.mat'));
save(SaveFullFileName, 'MAC_tracking_data_individual_pore_growth_videos_True_tracks');

% At this point, this data structure can either be loaded into:
% MACanalysis_TrackVideos_concatonate, to aggregate this data with other
% experiments, or, it can be loaded straight into the final script:
% MACanalysis_TrackVideos_tau_fitting, to calculate the time for individual
% pore growths, and plot histograms etc.











