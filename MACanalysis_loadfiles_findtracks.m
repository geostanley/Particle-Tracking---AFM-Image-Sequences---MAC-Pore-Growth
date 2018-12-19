%%%=== MACanalysis_loadfiles_findtracks ===%%%

% This is script 1/5.

% This script only works with the Nanoscope MATLAB utilities toolbox - 
% required to load the .spm image files.

% The purpose of this script is to load a video sequence, for the user to
% define a template (in this case a MAC pore), and then for the script to
% find things similar to the template throughout the sequence. If the object
% is found to be the same object from 1 image to the next, it creates 
% tracks of these things diffusing across the image. Whether a displaced
% object is indeed a track or a different object is defined by the tracking
% parameters (explained later).

% Therefore, when this script is run, it produces two videos. The first
% highlights objects found defined as being similar to the template, and
% the second shows the defined tracks. If either of these videos show
% incorrect results (false positives etc), either the threshold for the
% normalised cross-correlation, or the parameters defining the tracks,
% should be adjusted accordingly.

% The script runs as follows:

% It takes an image sequence specified by the user, applies a 1st order 
% plane fit to each image (setting the background to 0 nm at the same time),
% and displays the final image in the sequence to the user and asks for them
% to select a template. It then uses this template with the normxcorr2 function in the
% findfeatures_spmFiles function to find the coordinates of protrusions in an image.
% These coordinates are then sent into the simpletracker function to make
% 'tracks' of events captured throughout the video. This is ultimately
% stored in the tracking table: an Nx4 array with columns: x-coordinate,
% y-coordinate, frame numb, and track ID. This, along with the flattened
% image sequence (and some other parameters), are saved in a data
% structure, to be loaded back into the next script:
% MACanalysis_loadtracks_createTrackVideos, where videos of individual pore
% formation events are created.

%% Input data directory, file name, and output directory

clear variables
close all
clc

% Load data from DatFolder directory
DataFolder   = 'Z:\Users\George\Documents\PhD\MAC_Manuscript\MAC_Analysis_Everything\EPFL2_experiment\Data';
% Save outputs to:
outdirectory = 'Z:\Users\George\Documents\PhD\MAC_Manuscript\MAC_Analysis_Everything\EPFL2_experiment\DataStructures';

DataStructure_SaveName = 'MAC_EPFL2_tracktable';

% To review accuracy of tracking, set to 1. If incorrect, adjust tracking
% parameters and re-start script.
Show_video_with_tracks = 1; % set to 1 if would like to see

% File names and numbers
GenericFileName = 'c9';
File_Nos        = 0:465;

%%%=== Input imaging conditions ===%%%
NS_nine_one   = 'yes'; % if NS 9.1 or lower = 'yes'; if 9.2 or higher = 'no'
channel       = 1; % input channel number
ScanSize_nm   = 600; % input scan size in nm
scaleFactor_z = 0.2; % for getting correct z values (nm) (if microscope correctly calibrated, set to 1)
LineRate_Hz   = 39.1; % for calculating the frame rate later


% colour limits (z in nm) for images
clims = [0 15];

%% Tracking parameters: adjust to best find all pores

% this is the parameter for the normalised cross-correlation function.
threshold = 0.6; % Threshold for cross-correlation peaks (-1 to 1)
% if can see from video that pore counting is not good, run script again
% with adjusted threshold value

% Size of template image used in cross-correlation routine
CropSize_template_nm = 60; % (in nm)

%%%=== Tracking parameters ===%%%
% Maximum distance (in pixels) between two features in subsequent frames 
% to form a track
max_linking_distance = 30; % Something of the size of the template?
% Maximum number of frames over which a features be tracked if interrupted
% (e.g., if max_gap_closing = 3, a feature may be absent for 3-1=2 frames
% and still be considered as belonging to the same track
max_gap_closing = 5; % For longer movie with many frames, e.g. 5
% Minimum number of frames in which a feature needs to be detected to count
% as a track
mintracklength = 10; % A higher number, e.g. 10, helps to eliminate false counts, but requires a long sequence
% Show debugging info from simpletracker routine
debug = true;


%% Load data, re-calibrate the height, and store in cell array

% access the Nanoscope MATLAB toolbox for loading raw spm files 
NSMU = NSMatlabUtilities();

display('Creating file names and loading image sequence...')

% create cell array of filenames for image sequence
FileNames_Sequence = Create_FileNames_Cell_2(DataFolder, GenericFileName, File_Nos, NS_nine_one);
% pre-allocate cell array for loaded spm files (nm)
Vid_Seq_Mat = cell(size(File_Nos));
Vid_Seq_Raw = cell(size(File_Nos));

for i = 1:length(File_Nos)
    
    ImageFileName = FileNames_Sequence{i};
    NSMU.Open(ImageFileName);
    [heightdata, ~, ~] = NSMU.GetImageData(channel, NSMU.METRIC);
    
    % scale the height data to correct for microscope calibration file
    heightdata_scaled = heightdata .* scaleFactor_z;
    % save images into cell array
    Vid_Seq_Mat{i}    = heightdata_scaled;
    Vid_Seq_Raw{i}    = heightdata;
    
end

% take size of final image in series (this assumes all images in sequence
% are of the same size), and use this along with the LineRate_Hz to get the
% frame rate (ScanRate_Hz).
[row_img, col_img] = size(heightdata);
nmperpixel         = ScanSize_nm/col_img;
ScanRate_Hz        = row_img/LineRate_Hz;
    
%% Remove 1st order plane from each image in sequence (twice)
% this ensures all images are correctly processed and simultaneously sets
% the background SLB of each image to ~0 nm.

display('Apply 1st order plane fitting to each image in sequence...')

% pre-allocate cell array for saving flattened images
Vid_Seq_flat = cell(size(Vid_Seq_Mat));

for i = 1:length(Vid_Seq_Mat)
    
    % extract each image in sequence
    height_mat = Vid_Seq_Mat{i};

    % transform to XYZ array for plane fitting
    [XYZ_array_1] = Matrix_to_Nx3array(height_mat);

    % for plane fitting, take the bottom 50% of data
    BinWidth_nm      = 0.1; % the smaller the value the more accurate the XYZ indexing. 0.1nm should be fine.
    Plane_fit_mask_1 = 0.5;
    greater_than_1   = 0;
    
    % create new XYZ array of only bottom 50% of data
    [XYZ_array_for_plane_fit_1, ~] = XYZarray_indexed_by_percentage_height(XYZ_array_1, BinWidth_nm, Plane_fit_mask_1, greater_than_1);

    % find the plane of the indexed bottom 50% of height data
    [plane_1] = PlaneFit_XYZarray(height_mat, XYZ_array_for_plane_fit_1);
    % subtract plane from data. This operation also brings the
    % background = 0 nm.
    height_mat_1_flat = height_mat - plane_1;
        
    % now repeat this operation to ensure negligible background slope in
    % data
    [XYZ_array_2]    = Matrix_to_Nx3array(height_mat_1_flat);
    BinWidth_nm_2    = 0.1; % the smaller the value the more accurate the XYZ indexing. 0.1nm should be fine.
    Plane_fit_mask_2 = 0.3;
    greater_than_2   = 0;
    [XYZ_array_for_plane_fit_2, XYZ_array_remaining_data_2] = XYZarray_indexed_by_percentage_height(XYZ_array_2, BinWidth_nm_2, Plane_fit_mask_2, greater_than_2);
    [plane_2]         = PlaneFit_XYZarray(height_mat_1_flat, XYZ_array_for_plane_fit_2);
    height_mat_2_flat = height_mat_1_flat - plane_2;
    
    % save background flattened data into new cell array: height_data_flat
    Vid_Seq_flat{i} = height_mat_2_flat;
       
end

%% Select a template for the cross-correlation scores

[template] = SelectTemplate_even(Vid_Seq_flat{end}, CropSize_template_nm, nmperpixel, clims);

figure(), 
imagesc(template)
pbaspect([1 1 1]);
colormap(copper)
caxis([0 20])
title('Selected template for cross-correlation analysis', 'FontSize', 15)
xlabel('Pixels', 'FontSize', 13)
ylabel('Pixels', 'FontSize', 13)

%% Using normxcorr2, find number of completed MACs (for kinetics later)

%clear coordinates; clear imagedatalabelled; % Makes sure target cell arrays are empty
display('Analysing images...');

count = zeros(size(Vid_Seq_flat)); % Set count to zero
coordinates = cell(size(Vid_Seq_flat)); % pre-allocate coordinates cell


for i = 1:length(Vid_Seq_flat)
    
    % load each flattened image
    image = Vid_Seq_flat{i};
    
    % find shapes greater than threshold
    coordinates{i} = findfeatures_spmFiles(image, template, threshold);
    % save number of completed MACs
    count(i)       = length(coordinates{i});
        
end

%% Display video sequence showing when MACs are found

% if can see from video that pore counting is not good, run script again
% with adjusted threshold value (see top of script)

figure();
title('Number of completed MACs (count)', 'FontSize', 16)
pause on
hold on
for i = 1:length(Vid_Seq_flat)
    
    image = Vid_Seq_flat{i};
    
    
    % Each frame is plotted with the coordinate markers, for inspection
    imagesc(image);
    pbaspect([1 1 1]);
    colormap(copper);
    caxis([0 20]);
    title('Found protrusions', 'FontSize', 16)
    xlabel('Fast-scan axis (pix)')
    ylabel('Slow-scan axis (pix)')
    hold on;

    p = coordinates{i};
    if ~isempty(p) 
        plot(p(:,1),p(:,2),'ow','LineWidth',2,'MarkerSize',mean(size(template))+4);
    end;
    hold off;
    clear p;
    pause(0.02)
        
end

%% Plot number of completed MACs

figure(), plot(count, 'LineWidth', 2)
title('Number of completed MACs', 'FontSize', 16)
xlabel('Frame number', 'FontSize', 13)
ylabel('Count', 'FontSize', 13)

%% Next track features using simpletracker.m...

% First need to create newcoordinates set with all empty cells (i.e.,
% frames without features removed), to avoid error with tracking routines

nonemptyindices = find(~cellfun(@isempty,coordinates));% Find nonempty cells
newcoordinates  = coordinates(nonemptyindices); % Remove empty cells

firstframe = 1; % First image to be analysed - set to 1 if unknown
finalframe = length(newcoordinates); % Final image to be analysed

% Run tracking routine (see description in simpletracker.m)
[tracks, adjacency_tracks] = simpletracker(newcoordinates,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', debug);

% Put "empty frames" back again into track lists
% And remove tracks of features that appear
% in fewer than mintracklength frames
nanarray = NaN*ones(length(coordinates),1);
i=1;
while i<=length(tracks)
    testtrack=nanarray;
    testtrack(nonemptyindices)=tracks{i};
    if sum(~isnan(testtrack))<mintracklength
        tracks(i)=[]; % Delete track if shorter than mintracklength
    else
        tracks{i}=testtrack; % Insert NaN for frames not included in tracking
        i=i+1;
    end
end;

% Next put into matrix with 
% col1 x; col2 y; col3 framenumber; col4 trackID
trackingtable = createtracktable(coordinates,tracks);

% Recalculate track count (now too short tracks are eliminated)
trackingtable = createappearancetable(trackingtable,[0 length(Vid_Seq_flat)]);
for i=1:length(count)
    count(i)=length(find(trackingtable(:,3)==i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get the index values for each track ID and store in cell array
% This is to transform the large trackingtable into a more managable cell
% array later. The tracking table is an Nx4 array with columns: x, y, 
% frame numb, and track ID. This is reorganised into a cell array, where
% each cell is an array containing the index values for that track number.
% This can be used, if required, in the following script to get information
% on an entire track.

final_track_number = trackingtable(end,end);
track_numb_idx_cell = cell(1, final_track_number);

for i = 1:final_track_number
    track_numb_idx = find(trackingtable(:,4) == i);
    track_numb_idx_cell{i} = track_numb_idx;
    
end

%% Review accuracy of tracking

% Create images only covering analysed range

showimages = cell(1, length(firstframe:finalframe));
for i = firstframe:finalframe
    j = i - firstframe+1;
    showimages{j} = Vid_Seq_flat{i};
end;

% Plot results of tracking

if Show_video_with_tracks == 1

    % Determine number of tracks
    n_tracks = numel(tracks);
    colors   = linspecer(n_tracks); % create pretty colours for tracks

    % Make a list of all coordinates for all features
    % adjacency_tracks{i} contains the indices of track {i} in this list
    all_points = vertcat(coordinates{:});

    figure(5);
    hold on
    pause on

    for i = 1:length(showimages)

        figure(5);
        img = showimages{i};
        pbaspect([1 1 1]);
        imagesc(img);
        colormap(copper);
        caxis(clims);

        hold on

        for i_track = 1 : n_tracks

            % We use the adjacency tracks to retrieve the points coordinates. It
            % saves us a loop.

            track        = adjacency_tracks{i_track};
            track_points = all_points(track, :);
            plot(track_points(:,1), track_points(:, 2), 'Color', colors(i_track, :))

        end

        pause(0.005);
        hold off

    end

end

%% Save out important information for the plotting script

display('Saving built tracks to data structure...')

MAC_tracking_data.Vid_Seq_Raw          = Vid_Seq_Raw;
MAC_tracking_data.Vid_Seq_Original     = Vid_Seq_Mat;
MAC_tracking_data.showimages           = showimages;
MAC_tracking_data.CropSize_template_nm = CropSize_template_nm;
MAC_tracking_data.ScanSize_nm          = ScanSize_nm;
MAC_tracking_data.LineRate_Hz          = LineRate_Hz;
MAC_tracking_data.ScanRate_Hz          = ScanRate_Hz;
MAC_tracking_data.nmperpixel           = nmperpixel;
MAC_tracking_data.Vid_Seq_flat         = Vid_Seq_flat;
MAC_tracking_data.template             = template;
MAC_tracking_data.trackingtable        = trackingtable;
MAC_tracking_data.count                = count;
MAC_tracking_data.track_numb_idx_cell  = track_numb_idx_cell;
MAC_tracking_data.firstframe           = firstframe;
MAC_tracking_data.finalframe           = finalframe;
MAC_tracking_data.coordinates          = coordinates;
MAC_tracking_data.newcoordinates       = newcoordinates;
MAC_tracking_data.max_linking_distance = max_linking_distance;
MAC_tracking_data.max_gap_closing      = max_gap_closing;
MAC_tracking_data.mintracklength       = mintracklength;
MAC_tracking_data.tracks               = tracks;
MAC_tracking_data.adjacency_tracks     = adjacency_tracks;


SaveFullFileName = fullfile(outdirectory, strcat(DataStructure_SaveName, '.mat'));
save(SaveFullFileName, 'MAC_tracking_data');

% Re-load this data structure in the script: MACanalysis_loadtracks_createTrackVideos to
% make the library of videos of pore formation.