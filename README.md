# Particle tracking analysis of AFM image sequences reveals single-molecule kinetics of pore assembly by the membrane attack complex

The scripts (and corresponding functions) within this project are designed to analyse high-speed AFM image  sequences (.spm files). More specifically, they are designed to find features within an image sequence that resemble a template image, to build tracks of these features as they diffuse across the image sequences with time, and to crop out individual image sequences of these diffusing features, using the coordinates from the found tracks. 

In this study, these scripts are used to count the number of membrane attack complex (MAC) pores on a supported lipid bilayer within each image, to build tracks of these MAC pores as they grow and diffuse across the bilayer, to crop an image sequence of each found growing pore, and finally, to calculate how long it takes for each individual pore to grow.

The results from this analysis can be found here: https://www.biorxiv.org/content/early/2018/11/16/472274

There are five main scripts:

1. MACanalysis_loadfiles_findtracks.m (or MACanalysis_loadfiles_MATds_findtracks.m, which works without installing the NSMatlabUtilities toolbox) 
2. MACanalysis_loadtracks_createTrackVideos.m
3. MACanalysis_eliminateFalsetracks.m
4. MACanalysis_TrackVideos_concatonate.m
5. MACanalysis_TrackVideos_tau_fitting.m

Data sets to be used with the scripts can be provided upon request. 

The 5 scripts run as follows:

## Script 1: MACanalysis_loadfiles_findtracks.m

**MACanalysis_loadfiles_findtracks.m:** This script loads raw .spm image files (from Bruker AFM systems). (It therefore requires the installation of the Nanoscope MATLAB toolbox: NSMatlabUtilities. For this reason, another version of this script is provided that loads images directly from a MATLAB data structure: MACanalysis_loadfiles_MATds_findtracks.m. Read below for more details.) After loading the images, it applies a 1st order plane background subtraction to each image in the sequence. It shows the final image in the sequence (in which it is assumed complete MAC pores can be seen), and asks the user to select a template pore (this is done by clicking on the centre of one MAC and pressing enter). This template is used in a 2D cross-correlation routine to find similar features. The coordinates from these found feautures are entered into the simpletracker.m function to build tracks of particles diffusing across the image (the tracks are defined by the max_linking_distance and max_gap_closing parameters, defined at the top of the script). Tracks shorter than the mintracklength are then removed (this parameter is also set at the top of the script). All the tracks are stored in an Nx4 array (trackingtable), with columns: x-coordinate, y-coordinate, frame number, and track ID. This, along with the flattened image sequence (and some other parameters), are saved in a data structure, to be loaded back into the next script.

This script shows two video outputs. The first circles found protrusions that are similar to the template; the second shows the tracks. If these video outputs are not showing the correct results, the parameters must be adjusted and the script run again. These parameters include: the template size and the threshold value (for the 2D cross-correlation) used to find similar features; and the max_linking_distance, max_gap_closing, and mintracklength parameters used for building the tracks. 

Inputs: .spm files, imaging parameters (image size and line rate etc - see preamble).

Outputs: MATLAB data structure with flattened image sequence, imaging parameters, tracking table, and tracking parameters.

NB. This script only works with the Nanoscope MATLAB Utilities toolbox correctly installed. Another script is provided (MACanalysis_loadfiles_MATds_findtracks.m) which loads raw AFM images that have already been saved into  a MATLAB data structure (MAC_EPFL1_ImagesOnly_DataStructure_Files_0to478.mat). This allows for testing of the scripts without having to install the Nanoscope MATLAB toolbox.

## Script 2: MACanalysis_loadtracks_createTrackVideos.m

**MACanalysis_loadtracks_createTrackVideos.m:** This script loads the data structure output from the previous script. It then reorganises the tracking table: for a given track, it takes the coordinates of when the feature is first seen to appear and adds in frames going backwards in time (the number of frames is defined by preframe); it then also takes frames going forward in time, using the coordinates of the track (the number of frames is defined by postframe). This enables us to, later, crop out image sequences of growing pores, from background lipid bilayer to fully grown MAC pores. The reorganised tracking table (defined as trackappearance) is then used to produce an MP4 movie file of the tracked MACs throughout the image sequence. The trackappearance table is turned into a cell array in which each array contains the tracking information of a growing pore. This is used to make a new cell array (Cropped_pore_growth_cell), containing the cropped image sequences of each track (the size of these images is the same size as the template image, defined in the previous script). The mean height of each frame within the cropped image sequences is also calculated (this is later used as a measure of pore completion). MP4 files are also made of these cropped tracks. These are useful later for reviewing the accuracy of the tracking (and ultimately, for eliminating false positives etc).

At this point, have made a library of cropped image sequences from the tracks (using the preframe and postframe parameters and subsequent reorganisation of the track table), and have calcualted the mean height profile for each cropped image sequence. This information is saved into another MATLAB data structure.

Input: MATLAB data structure produced from 'MACanalysis_loadfiles_findtracks.m'.

Outputs: MATLAB data structure containing the cropped image sequences (Cropped_pore_growth_cell), the mean height arrays for each sequence, the template size, and the imaging frame rate (for producing time vectors etc).

## Script 3: MACanalysis_eliminateFalsetracks.m

**MACanalysis_eliminateFalsetracks.m:** It must be remembered that the ultimate aim of these scripts is to calculate the time taken for individual MAC pores to grow, and that the tracking of these pores is done towards this aim. However, tracks cannot distinguish between complete MACs that have drifted into the field-of-view during imaging, and MAC pores that have oligomerised within the field-of-view during the experiment. This is the reason for saving out all the MP4 files produced in the previous script. These must now be reviewed, and  the track numbers corresponding to real pore growth events noted down. These numbers must be entered into the array "True_tracks_idx" within script 3 (MACanalysis_eliminateFalsetracks.m). The script will then load the MATLAB data structure from the previous script, but only save out the information for the track numbers specified in "True_tracks_idx". These are again saved out into a new data structure.

Inputs: MATLAB data structure produced from 'MACanalysis_loadtracks_createTrackVideos.m', and desired track numbers (entered into the "True_tracks_idx" array).

Output: MATLAB data structure.

## Script 4: MACanalysis_TrackVideos_concatonate.m

**MACanalysis_TrackVideos_concatonate.m:** Before calculating the time for each pore growth event (and plotting histograms etc), it is possible to concatonate the results from several independent experiments (if desired - this step is not necessary and can be skipped). Simply input the directories and data structure filenames of the experiments to be aggregated.

Inputs: Directories and filenames of data structures to be concatonated.

Output: Concatonated data structure.

## Script 5: MACanalysis_TrackVideos_tau_fitting.m

**MACanalysis_TrackVideos_tau_fitting.m:** After elinating false tracks (and perhaps concatonating different experiments), the mean height cell vectors from each growing pore sequence are used to calculate time for pore growth. The mean height vectors are filtered (to account for image noise) and normalised between 0 and 1 (arbitrary units). Each normalised height vector is then fit with a sinusoidal fitting equation to calculate the width of the transition from background lipid bilayer to saturated mean height (pore completion). If the width of this transition (and hence time of reaction coordinate) is found to be longer than the image sequence, or t_0 is negative, these are considered bad fits and the data is removed. Histograms of time for individual pores to form are then produced, along with videos of their pore growth and concomitant mean height profiles. The mean time for individual pore formation is also reported.

NB: This script can either load a data structure straight from 3 (MACanalysis_TrackVideos_tau_fitting.m), or, from 4 (MACanalysis_TrackVideos_concatonate.m). The user must input whether this is a concatonated data set or not.(Concatonated_data_set=0 or 1).

Data structures are provided which can be directly loaded into this final script for testing: "MAC_EPFL1_pore_track_vids_TrueTracks.mat", which hold the results from the data in "MAC_EPFL1_ImagesOnly_DataStructure_Files_0to478.mat"; and "MAC_EPFL1_EPFL2_concatonated.mat", which are the aggregated results from two independent experiments.
