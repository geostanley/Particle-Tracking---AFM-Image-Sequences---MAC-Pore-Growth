%%%=== MACanalysis_loadtracks_createTrackVideos ===%%%

% This is script 2/5.

% In the previous script, raw images were loaded and the images were
% flattened. A template of a MAC pore was chosen from the final frame in
% the video sequence, and then the MAC pores were tracked using a 2D
% cross-correlation routine. The tracking results were saved out as a 
% MATLAB data structure.

% In this script, the data structure is loaded, and a new tracking table is
% created. This table will be used to create videos of the growth of each
% pore individually. However, it further incorporates frames from the image
% sequence both from previous to the first frame the MAC pore was found 
% (using those coordinates); and after the pore was completed. This enables
% a later analysis of the pore growth from baseline SLB, to saturation of
% pore completion.

% After the completion of the new tracking table, each track is organised
% into a new cell array, in which each array is the cropped image sequence
% of an individual pore growing.

% The mean height of each frame of the individual pore growth videos are
% calculated (these are later used to calculate time for pore completion).
% All data is again saved out as a MATLAB data structure.

% The count data is also plotted against time (s) (and the plot saved).

%% Input data directory, file name, and output directory

clear variables
close all
clc

Load_directory          = 'Z:\Users\George\Documents\PhD\MAC_Manuscript\MAC_Analysis_Everything\EPFL2_experiment\DataStructures';
outdirectory            = 'Z:\Users\George\Documents\PhD\MAC_Manuscript\MAC_Analysis_Everything\EPFL2_experiment\DataStructures';
outdirectory_vids_plots = 'Z:\Users\George\Documents\PhD\MAC_Manuscript\MAC_Analysis_Everything\EPFL2_experiment\FiguresVideos';

DataStructure_LoadName = 'MAC_EPFL2_tracktable';
DataStructure_SaveName = 'MAC_EPFL2_pore_track_vids';

%%%=== Choose save names for output of videos and plots ===%%%
Video_Seq_title   = 'MAC_ecolilipid_32kHzPFT_EPFL2'; % for entire, annotated video sequence
countfilename     = 'MAC_ecolilipid_32kHzPFT_EPFL2-count.txt'; % for count plot
countfilename_fit = 'MAC_ecolilipid_32kHzPFT_EPFL2-count_fit.txt'; % for count plot with fitting
parsfilename      = 'MAC_ecolilipid_32kHzPFT_EPFL2-pars.txt'; % Set name for parameter file

%%%=== For viewing and saving out the movies ===%%%
Save_whole_movie_seq           = 1;
Save_individual_pore_movie_seq = 1;
plane_fitting_cropped_img      = 1; % 1st order plane subtraction to each cropped frame


%%%=== Parameters for the mp4 video appearance ===%%%
% specify how many frames/sec you wish the mp4 files to be saved at
Saved_frames_per_sec = 6; 
% enter colour limits (nm) for appearance of videos
clims = [0 15];
% size of template in nm
% CropSize_template_nm = 60;
Radial_crop_nm = 25;

%%%=== Parameters for length of tracks pre- and post-pore formation ===%%%
% these inputs are for the appearance of the individual pore videos and
% also for their kinetic analysis of formation. It is the number of frames
% (in s) before and after the first track appeared. You probably want the
% postfram_s to be longer than pre, as you must account for the time taken
% to form the pore, as well as wanting some information after the pore's
% formation.
preframe_s  = 200;
postframe_s = 550;

% nice colours for plotting later
colours = linspecer(4);
colour1 = colours(1, :);
colour2 = colours(2, :);
colour3 = colours(3, :);
colour4 = colours(4, :);

%% Load the data structure

display('Loading tracking data structure...')
LoadFullFileName = fullfile(Load_directory, DataStructure_LoadName);
load(strcat(LoadFullFileName, '.mat'));

Vid_Seq_Original        = MAC_tracking_data.Vid_Seq_Original;
Vid_Seq_flat            = MAC_tracking_data.Vid_Seq_flat;
showimages              = MAC_tracking_data.showimages;
trackingtable           = MAC_tracking_data.trackingtable;
count                   = MAC_tracking_data.count;
track_numb_idx_cell     = MAC_tracking_data.track_numb_idx_cell;
template                = MAC_tracking_data.template;
FrameRate_s             = MAC_tracking_data.ScanRate_Hz;
firstframe              = MAC_tracking_data.firstframe;
finalframe              = MAC_tracking_data.finalframe;
coordinates             = MAC_tracking_data.coordinates;
newcoordinates          = MAC_tracking_data.newcoordinates;
max_linking_distance    = MAC_tracking_data.max_linking_distance;
max_gap_closing         = MAC_tracking_data.max_gap_closing;
mintracklength          = MAC_tracking_data.mintracklength;
adjacency_tracks        = MAC_tracking_data.adjacency_tracks;
tracks                  = MAC_tracking_data.tracks;
CropSize_template_nm    = MAC_tracking_data.CropSize_template_nm;
ScanSize_nm             = MAC_tracking_data.ScanSize_nm;

%% Set names for movie file of entire video sequence

% change preframe and postframe from secs back to frame numbers
preframe  = round(preframe_s/FrameRate_s); % Number of frames over which to mark position of new track (prior to its appearance)
postframe = round(postframe_s/FrameRate_s); % Number of frames over which to mark newly appeared track

% Create matrix in tracking entry format, but only selecting the first
% "postframe" entries for each track, and inserting "preframe" entries with
% same coordinates as first entry. This is to facilitate highlighting the
% appearance of new tracks

% In brief, this reorganises the tracking table. For a given track, it
% starts with the coordinates of when it is first seen to appear. Using
% these coordinates, it then adds in frames (for preframe number of
% frames), going backwards in time. This should be background SLB. From
% this point (of first appearance), it also takes postframe number of frames 
% forward (using the coordinates of the track). This enables us to, later, 
% crop out image sequences of growing pores, from background SLB to fully 
% grown MAC pores. If, for a given track, preframe or postframe exceeds
% either the start or end of the image sequence, the start or end of the
% sequence is used.
trackappearance = createappearancetable(trackingtable,[preframe postframe]);
%
% Export table with parameter settings for analysis
pars=table([ firstframe; finalframe; max_linking_distance; max_gap_closing; mintracklength; preframe; postframe ], 'RowNames',{'firstframe','finalframe','max_linking_distance','max_gap_closing','mintracklength','preframe','postframe'});
writetable(pars,fullfile(outdirectory,parsfilename),'Delimiter',' ','WriteRowNames',true);

% Export table with count data
dlmwrite(fullfile(outdirectory,countfilename),count);

%% Make movie of whole sequence with numbered tracked events

% Plot tracking results with numbered tracks

if Save_whole_movie_seq == 1;

    Vid_entire_Seq_Name = strcat(Video_Seq_title, '_Entire_Tracking_Vid_Seq');
    Video_entire_seq_Save_FullFileName = fullfile(outdirectory_vids_plots, Vid_entire_Seq_Name);

    vid_entire = VideoWriter(Video_entire_seq_Save_FullFileName, 'MPEG-4');
    vid_entire.FrameRate = 6;

    figure();
    hold on;
    open(vid_entire);

    for i=1:length(showimages)
        
        entire_img_flat = showimages{i};  
        imagesc(entire_img_flat, clims); 
        colormap(copper);
        caxis(clims)
        pbaspect([1 1 1]);
        hold on; % Plot image in figure
        coordinateindices = find(trackappearance(:,3) == i);
        if ~isempty(coordinateindices)
            imagecoordinates = trackappearance(coordinateindices,:);
            
            for j=1:length(coordinateindices)
                p=imagecoordinates(j,1:2);
                ptext=p+0.75*size(template(:,:,1));
                n=imagecoordinates(j,4);
                plot(p(:,1),p(:,2),'sw','LineWidth',2,'MarkerSize',1.5*mean(size(template)));
                text(ptext(1),ptext(2),num2str(n),'FontSize',16,'Color',colour2);
            end;
            
        end
        set(gca, 'xtick', [])
        set(gca, 'ytick', [])
        hold off;
        h=getframe; % Convert figure into image
        writeVideo(vid_entire, h)
        
    end
    
    close(vid_entire);
    
end

%% Plot count of features as function of framenumber

time_s = linspace(FrameRate_s, FrameRate_s*length(count), length(count));

figure();
plot(time_s, count, 's','MarkerEdgeColor', colour3, 'MarkerFaceColor', colour3, 'MarkerSize', 4, 'Linewidth', 0.2)
ylim([0 50])
xlim([0 3500])
title('Number of pores as a function of time')
set(gca, 'FontSize', 13)
xlabel('Time (s)')
ylabel('Number of pores')

%% Rate of MAC appearance - fit A x (1 - exp(-t/tau)) to counts

% copy time and count vectors for fitting
time_s_C9               = time_s;
height_data_for_fitting = count;

%%%=== Fit the equation

% suppress the 'Optimization complete' message
opts = optimset('Display', 'off');

% initial conditions for b(1) and b(2) respectively
ic = [mean(height_data_for_fitting) 1000];

%%%=== Form of the fitting equation ===%%% 
%%%=== A x (1 - exp(-t/tau))        ===%%%
Kinetics_Fit_Equation = @(b, time_s_C9) ( b(1)* (1 - exp(-time_s_C9/b(2))) );

% Use a non-linear least squares regression function
[C,R,J,CovB,MSE] = nlinfit(time_s_C9, height_data_for_fitting, Kinetics_Fit_Equation, ic, opts);

% save out the fitted parameters
A   = C(1);
tau = C(2);

% get the 95% confidence interval
[ypred, delta] = nlpredci(Kinetics_Fit_Equation, time_s_C9, C, R,'Covar', CovB,...
                         'MSE', MSE, 'SimOpt', 'on');
lower = ypred - delta;
upper = ypred + delta;

xax = time_s_C9;
figure(); 
plot(time_s_C9, height_data_for_fitting, 's','MarkerEdgeColor', colour3, 'MarkerFaceColor', colour3, 'MarkerSize', 4, 'Linewidth', 0.2)
ylim([0 50])
xlim([0 3500])
set(gca, 'LineWidth', 2, 'FontSize', 13)
hold on 
plot(xax, Kinetics_Fit_Equation(C,xax), 'Color', colour2, 'Linewidth', 2) 
hold on
plot(time_s_C9, lower, 'LineStyle', '--', 'Color', colour2)
hold on
plot(time_s_C9, upper, 'LineStyle', '--', 'Color', colour2)
ylabel('Count', 'FontSize', 14)
xlabel('Time (s)', 'FontSize', 14)
title(['A = ', num2str(A), '; \tau = ', num2str(tau)], 'FontSize', 16)

fig_name = strcat('MAC_count', '_tau_fit_confidence_interval');
saveas(gca, fullfile(outdirectory_vids_plots, fig_name), 'pdf')


%% Plotting tracks from track appearance table

% transform the trackappearance table into a cell-array (for ease of manipulation later)
final_track_number = trackappearance(end,end);
track_app_table_cell = cell(1, final_track_number);

for i = 1:final_track_number
    
    track_numb_app_idx = find(trackappearance(:,4) == i);  
    track_app_table_cell{i} = trackappearance(track_numb_app_idx, :);
    
end

%% Crop individual tracked pore formation events and save into cell array

% create a circular mask to remove neighbouring pores in individual pore
% growth image sequences. This is for calculating the mean height of the
% frames later.
[~, CircleMatrix, ~]         = RadialBins(CropSize_template_nm, template, 5, 5, 14);
CircleMaskLogic              = CircleMatrix <= Radial_crop_nm;
CircleMaskLogic_pixel_number = sum(CircleMaskLogic(:));
template_masked              = template .* CircleMaskLogic;

% show example size of circular mask on the template used for the 2D
% normalised cross-correlation
figure();
imagesc(template_masked)
colormap('copper')
caxis(clims)
title('Circle mask of template (example)', 'FontSize', 16)
xlabel('Pixels', 'FontSize', 13)
ylabel('Pixels', 'FontSize', 13)
pbaspect([1 1 1])
set(gca, 'FontSize', 13, 'LineWidth', 2)

%% Create the cell array of individual pore growth videos and calculate mean heights

% take size of template (used for normxcorr2) for padding (see later)
[row_temp, col_temp] = size(template);
% create a Crop_Size_Half variable for imcrop function
Crop_Size_Half = round(col_temp/2);
% pre-allocate cell array for image sequences of cropped pores
Cropped_pore_growth_cell = cell(size(track_app_table_cell));
mean_height_cell         = cell(size(track_app_table_cell));

for n=1:length(Cropped_pore_growth_cell)
    
    % pull-out x,y coordinates, and frame numbers for individually tracked
    % pores
    coordinate_frame_numb_seq = track_app_table_cell{n};
    % create cell-array of same size to store cropped images
    cropped_img_seq = cell(1, length(coordinate_frame_numb_seq(:,end)));
    
    mean_height_array = zeros(length(coordinate_frame_numb_seq), 2);
    
    for i=1:length(coordinate_frame_numb_seq)
     
        % x,y coordiantes plus size of padding
        centre_x = coordinate_frame_numb_seq(i,1)+row_temp;
        centre_y = coordinate_frame_numb_seq(i,2)+col_temp;
               
        % go to top left for imcrop
        crop_x = centre_x - Crop_Size_Half;
        crop_y = centre_y - Crop_Size_Half;
        
        % pull-out correct frame number and load entire image
        entire_img = Vid_Seq_flat{coordinate_frame_numb_seq(i,3)};
        [row_img, col_img] = size(entire_img);

        % pad image with zeroes around edge for cropping of events at edge
        padded_img = zeros(row_img+(2*row_temp), col_img+(2*col_temp));
        padded_img(row_temp+1:end-row_temp, col_temp+1:end-col_temp) = entire_img;
        
        % crop
        height_mat = imcrop(padded_img, [crop_x crop_y (row_temp-1) (col_temp-1)]);
        
        % 1st order plane background subtraction to each image
        if plane_fitting_cropped_img == 1
            % 1st order plane fit cropped images to ensure flat and background
            % at 0nm
            % transform to XYZ array for plane fitting
            [XYZ_array_1] = Matrix_to_Nx3array(height_mat);

            % find 1sr order plane of entire image
            BinWidth_nm = 0.1;
            Plane_fit_mask_1 = 1;
            greater_than_1   = 0;
            % index XYZ array
            [XYZ_array_for_plane_fit_1, ~] = XYZarray_indexed_by_percentage_height(XYZ_array_1, BinWidth_nm, Plane_fit_mask_1, greater_than_1);
            % get plane and subtract from height data
            [plane_1] = PlaneFit_XYZarray(height_mat, XYZ_array_for_plane_fit_1);
            height_mat_2 = height_mat - plane_1;

            % now repeat but only apply plane fitting to bottom 50% of data (background)
            [XYZ_array_2] = Matrix_to_Nx3array(height_mat_2);

            Plane_fit_mask_1 = 0.5;
            greater_than_1   = 0;

            [XYZ_array_for_plane_fit_2, ~] = XYZarray_indexed_by_percentage_height(XYZ_array_2, BinWidth_nm, Plane_fit_mask_1, greater_than_1);
            [plane_2] = PlaneFit_XYZarray(height_mat_2, XYZ_array_for_plane_fit_2);
            cropped_img_flat = height_mat_2 - plane_2;
        else
            cropped_img_flat = height_mat;
        end
        
        % crop out a circle around the pore
        cropp_img_flat_masked = cropped_img_flat .* CircleMaskLogic;
        
        % take mean height
        mean_height = sum(cropp_img_flat_masked(:))/CircleMaskLogic_pixel_number;
        mean_height_array(i, 2) = mean_height;   
        
        % save into cell-array
        cropped_img_seq{i} = cropp_img_flat_masked;

    end
    
    % create time vector for mean height array
    mean_height_array(:, 1) = linspace(0, FrameRate_s*length(mean_height_array), length(mean_height_array));
    
    % save into cell array
    mean_height_cell{n} = mean_height_array;    
    % save image sequence
    Cropped_pore_growth_cell{n} = cropped_img_seq;
    
end

%% Create movie sequences of individual pore formation events - and save

if Save_individual_pore_movie_seq == 1

    for n = 1:length(Cropped_pore_growth_cell)

        Cropped_img_seq = Cropped_pore_growth_cell{n};

        Seq_numb = num2str(n);
        Video_Seq_title_savename = strcat(Video_Seq_title, '_', Seq_numb);
        Video_Save_FullFileName = fullfile(outdirectory_vids_plots, Video_Seq_title_savename);

        vid = VideoWriter(Video_Save_FullFileName, 'MPEG-4');
        vid.FrameRate = Saved_frames_per_sec;
        figure();
        hold on
        open(vid);
        for i=1:length(Cropped_img_seq)

            imagesc(Cropped_img_seq{i}, clims);
            colormap('copper');
            pbaspect([1 1 1]);
            set(gca, 'xtick', [])
            set(gca, 'ytick', [])
            Img_frame = getframe;
            writeVideo(vid, Img_frame);

        end
        hold off
        close(vid);
        close(gcf);

    end

end

%% Plot mean heights of image sequences

figure(),
hold on
for i = 1:length(mean_height_cell)

    
    mean_height_array = mean_height_cell{i};
    time_vector       = mean_height_array(:,1);
    ave_height_vector = mean_height_array(:,2);
    plot(time_vector, ave_height_vector, 'x', 'Color', colour3)
    
end
title('Individual pore growths (t_{0} not yet properly defined)', 'FontSize', 16)
xlabel('Time (s)', 'FontSize', 13)
ylabel('Mean frame height (nm)', 'FontSize', 13)
set(gca, 'LineWidth', 2, 'FontSize', 13)
axis([0 800 0 6])

%% turn tracking table into cell using track_numb_idx_cell

trackingtable_cell = cell(1, length(track_numb_idx_cell));
for i = 1:length(track_numb_idx_cell)
    
    track_idx             = track_numb_idx_cell{i};
    track_info            = trackingtable(track_idx, :);
    trackingtable_cell{i} = track_info;
    
end

%% Save out individual track videos with mean height vectors

% This is so the data can later be collated with other experiments to
% produce histograms of individual pore growth times from several imaging
% experiments

% First, however, each growing pore must be reviewed manually and
% cross-checked with the video of all pores being found. This is to check a
% video of an individual pore is indeed just one pore, and the tracking has
% not moved from one pore to another. This is done in the following script:
% MACanalysis_eliminateFalsetracks

display('Saving individual pore growth sequences into data structure...')

MAC_tracking_data_individual_pore_growth_videos.FrameRate_s              = FrameRate_s;
MAC_tracking_data_individual_pore_growth_videos.CropSize_template_nm     = CropSize_template_nm;
MAC_tracking_data_individual_pore_growth_videos.mean_height_cell         = mean_height_cell;
MAC_tracking_data_individual_pore_growth_videos.Cropped_pore_growth_cell = Cropped_pore_growth_cell;
MAC_tracking_data_individual_pore_growth_videos.trackingtable_cell       = trackingtable_cell;


SaveFullFileName = fullfile(outdirectory, strcat(DataStructure_SaveName, '.mat'));
save(SaveFullFileName, 'MAC_tracking_data_individual_pore_growth_videos');

% Now have created library of individual pore growth videos. Review the MP4
% files to check for false positives etc that you do not wish to carry
% forward into the analysis. Load this data structure into:
% MACanalysis_eliminateFalsetracks to eliminate these false positives.