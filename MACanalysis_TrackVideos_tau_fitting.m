%%%=== MACanalysis_TrackVideos_tau_fitting ===%%%

% This is script 5/5.

% This script loads either the concatonated data structure, or a data
% structure straight from |MACanalysis_eliminateFalsetracks|. It takes the
% mean frame height from each pore growth video, filters it, normalises
% it, and then fits it with a geometric equation (defined later). This is
% used to define the width of the transition from background SLB to
% completed pore, and hence the time for individual pore growth. Videos are
% created with these plots, and a histogram of pore growth times are made
% and saved out.

%% Input data directory, file name, and output directory

clear variables
close all
clc

%%%=== Enter load directory and file name ===%%%
Load_directory         = 'Z:\Users\George\Documents\PhD\MAC_Manuscript\MAC_Analysis_Everything\Concatonated_data_sets\DataStructures';
DataStructure_LoadName = 'MAC_EPFL1_EPFL2_concatonated';

% Is this data structure being loaded from the concatonate script (==1), or
% straight from |MACanalysis_eliminateFalsetracks| (==0)?
Concatonated_data_set = 1;

%%%=== Enter output directory for videos and plots ===%%%
outdirectory       = 'Z:\Users\George\Documents\PhD\MAC_Manuscript\MAC_Analysis_Everything\Concatonated_data_sets\FiguresVideos';
outdirectoryimgseq = 'Z:\Users\George\Documents\PhD\MAC_Manuscript\MAC_Analysis_Everything\Concatonated_data_sets\IndividualPoreImgSeq';
Video_Seq_title_pg = 'EPFL1_EPFL2_cat';

%%%=== Do you want to save out plots and plotting videos? ===%%%
save_tau_fits_plots  = 1;
save_plotting_videos = 1;
Saved_frames_per_sec = 6; % for frame rate of saved videos (~6 recommended)
save_pore_as_image   = 1; % if want to save image sequence of a particular pore for figures (0 or 1)
pore_n               = 5; % select desired pore for saving all images in sequence

% histogram bin widths (s)
binwidth_s = 60;

% colour limits for individual pore growth videos
clims = [0 15];

% (for plotting with nice colours later)
C = linspecer(4);
c1 = C(1, :);
c2 = C(2, :);
c3 = C(3, :);
c4 = C(4, :);

%% Load the data structure

display('Loading tracking data structure...')
LoadFullFileName = fullfile(Load_directory, DataStructure_LoadName);
load(strcat(LoadFullFileName, '.mat'));

if Concatonated_data_set == 1 % if loaded from concatonated data struture

    img_seq              = MAC_True_tracks_concatonated.Cropped_pore_growth_cell_true_tracks;
    mh_track_all         = MAC_True_tracks_concatonated.mean_height_true_tracks_cell;
    firsttrack_sec_array = MAC_True_tracks_concatonated.firsttrack_sec_array;
    template_nm_array    = MAC_True_tracks_concatonated.template_nm_array;
    
else % if loaded straight from |MACanalysis_eliminateFalsetracks|
    
    img_seq              = MAC_tracking_data_individual_pore_growth_videos_True_tracks.Cropped_pore_growth_cell_true_tracks;
    mh_track_all         = MAC_tracking_data_individual_pore_growth_videos_True_tracks.mean_height_true_tracks_cell;
    firsttrack_sec_array = MAC_tracking_data_individual_pore_growth_videos_True_tracks.firsttrack_sec_array;
    CropSize_template_nm = MAC_tracking_data_individual_pore_growth_videos_True_tracks.CropSize_template_nm;
    
    % create an array of template sizes for plotting in the final script
    template_array = ones(size(img_seq));
    template_nm_array = template_array .* CropSize_template_nm;
    
end

%% Filter data mean height arrays
% this is to reduce interference from imaging noise in the tau fitting
% later

mh_track_sg  = cell(size(mh_track_all));

for i = 1:length(mh_track_all)
    
    t_vec = mh_track_all{i}(:,1);
    h_vec = mh_track_all{i}(:,2);
    
    h_vec_sg = smooth(h_vec, 'sgolay', 1);   
    mh_track_sg{i} = [t_vec, h_vec_sg];
    
end

%% Normalise the filtered, mean height arrays

mh_track_sg_norm  = cell(size(mh_track_sg));

for i = 1:length(mh_track_sg_norm)
    
    t_vec_sg = mh_track_sg{i}(:,1);
    h_vec_sg = mh_track_sg{i}(:,2);
    
    % normalise by taking the final 5 values (assumed to be
    % saturation), averaging, and dividing by this value
    h_saturation  = sum(h_vec_sg(end-4:end))/5;
    h_vec_sg_norm = h_vec_sg ./ h_saturation;  
    
    mh_track_sg_norm{i} = [t_vec_sg, h_vec_sg_norm];  
    
end

%% For fitting, will use the filtered, normalised mean height arrays
% this is written in explicitly here such that it is easy to change and
% fit to the original arrays, or the filtered arrays before normalisation.

mh_track_fitting = mh_track_sg_norm;

%% Fit geometric equation

% save out values for tau, 2*tau, and 3*tau for comparison if wanted later

% pre-allocate
Params_cell = cell(1,length(mh_track_fitting));

delta_t_tau_sg  = zeros(size(mh_track_fitting));
delta_t_2tau_sg = zeros(size(mh_track_fitting));
delta_t_3tau_sg = zeros(size(mh_track_fitting));

t0_tau_sg  = zeros(size(mh_track_fitting));
t0_2tau_sg = zeros(size(mh_track_fitting));
t0_3tau_sg = zeros(size(mh_track_fitting));

tend_tau_sg  = zeros(size(mh_track_fitting));
tend_2tau_sg = zeros(size(mh_track_fitting));
tend_3tau_sg = zeros(size(mh_track_fitting));

% Form of the equation - simple geometric equation to define width of
% transition
Fitting_eqn = @(b,t) (b(1) * tanh((t - b(2))/ b(3)) + b(4));

% initial conditions
ic = [0.5, 100, 15, 2];

for i = 1:length(mh_track_fitting)

    % extract time and height vectors
    t        = mh_track_fitting{i}(:,1);
    h_vec_sg = mh_track_fitting{i}(:,2);
    
    % non-linear least-squares regression fit
    Params = nlinfit(t, h_vec_sg, Fitting_eqn, ic);
    
    % save out fit parameters
    Params_cell{i} = Params;
    
    % parameters
    A         = Params(1);
    t_centre  = Params(2);
    tau       = Params(3);
    B         = Params(4);
    
    % save tau, 2tau, 3tau, t_0 and t_end   
    tau_2     = 2*tau;
    tau_3     = 3*tau;
    
    t0   = t_centre - 1.5*tau;
    tend = t_centre + 1.5*tau;

    % for plotting of lines later
    t0_x   = [t0 t0];
    t0_y   = [0 round(max(h_vec_sg))+0.5];
    tend_x = [tend tend];
    tend_y = [0 round(max(h_vec_sg))+0.5];
            
    delta_t_tau_sg(i)  = tau;
    delta_t_2tau_sg(i) = tau_2;
    delta_t_3tau_sg(i) = tau_3;
    
    t0_tau_sg(i)  = t_centre - 0.5*tau;
    t0_2tau_sg(i) = t_centre - 0.5*tau_2;
    t0_3tau_sg(i) = t_centre - 0.5*tau_3;

    tend_tau_sg(i)  = t_centre + 0.5*tau;
    tend_2tau_sg(i) = t_centre + 0.5*tau_2;
    tend_3tau_sg(i) = t_centre + 0.5*tau_3;
    
end

%% After completing the fitting, find bad fits and remove.

% remove data if t_0 from fit is negative, or if calculated reaction time
% is greater than length of data

keep_track = zeros(size(delta_t_3tau_sg));

for i = 1:length(delta_t_3tau_sg)
    
    Dt_3tau = delta_t_3tau_sg(i);
    t0_3tau = t0_3tau_sg(i);
    
    max_time = max(mh_track_fitting{i}(:,1));
    
    if t0_3tau<0 || Dt_3tau>max_time      
        keep_track(i) = 0;       
    else      
        keep_track(i) = 1;      
    end
    
end

keep_tracks_idx = find(keep_track == 1);

%% Go through data and keep selected tracks

numb_selected_fits = sum(keep_track);

img_seq_final          = cell(1, numb_selected_fits);
mh_track_fitting_final = cell(1, numb_selected_fits);
Params_cell_final      = cell(1, numb_selected_fits);
template_nm_final      = zeros(1, numb_selected_fits);

delta_t_tau_sg_final  = zeros(1, numb_selected_fits);
delta_t_2tau_sg_final = zeros(1, numb_selected_fits);
delta_t_3tau_sg_final = zeros(1, numb_selected_fits);

t0_tau_sg_final  = zeros(1, numb_selected_fits);
t0_2tau_sg_final = zeros(1, numb_selected_fits);
t0_3tau_sg_final = zeros(1, numb_selected_fits);

tend_tau_sg_final  = zeros(1, numb_selected_fits);
tend_2tau_sg_final = zeros(1, numb_selected_fits);
tend_3tau_sg_final = zeros(1, numb_selected_fits);

firsttrack_sec_array_final = zeros(1, numb_selected_fits);

for i = 1:length(keep_tracks_idx)
    
    track_numb = keep_tracks_idx(i);
            
    img_seq_final{i}          = img_seq{track_numb};
    mh_track_fitting_final{i} = mh_track_fitting{track_numb};
    Params_cell_final{i}      = Params_cell{track_numb};
    template_nm_final(i)      = template_nm_array(track_numb);

    delta_t_tau_sg_final(i)  = delta_t_tau_sg(track_numb);
    delta_t_2tau_sg_final(i) = delta_t_2tau_sg(track_numb);
    delta_t_3tau_sg_final(i) = delta_t_3tau_sg(track_numb);

    t0_tau_sg_final(i)  = t0_tau_sg(track_numb);
    t0_2tau_sg_final(i) = t0_2tau_sg(track_numb);
    t0_3tau_sg_final(i) = t0_3tau_sg(track_numb);

    tend_tau_sg_final(i)  = tend_tau_sg(track_numb);
    tend_2tau_sg_final(i) = tend_2tau_sg(track_numb);
    tend_3tau_sg_final(i) = tend_3tau_sg(track_numb);
    
    firsttrack_sec_array_final(i) = firsttrack_sec_array(track_numb);
    
end

%% Plot tau fits for accepted tracks

% there is some repetition and redundancy in this for loop but it's just
% for plotting and it works just fine so I can't be bothered to improve it.
% So there we are!
    
for i = 1:length(delta_t_3tau_sg_final)
    
    % extract time and height vectors
    t        = mh_track_fitting_final{i}(:,1);
    h_vec_sg = mh_track_fitting_final{i}(:,2);
    
    Params = Params_cell_final{i};
    
    t_centre  = Params(2);
    tau       = Params(3);
    
    tau_2     = 2*tau;
    tau_3     = 3*tau;
    
    t0   = t_centre - 1.5*tau;
    tend = t_centre + 1.5*tau;

    % for plotting of lines
    t0_x   = [t0 t0];
    t0_y   = [0 round(max(h_vec_sg))+0.5];
    tend_x = [tend tend];
    tend_y = [0 round(max(h_vec_sg))+0.5];
    
    t0_x_plot = linspace(t0, t0, 200);
    t0_y_plot = linspace(0, 1.5, 200);
    tend_x_plot = linspace(tend, tend, 200);
    tend_y_plot = linspace(0, 1.5, 200);

    ax1 = figure();
    plot(t, h_vec_sg, 's', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5])
    hold on
    plot(t, Fitting_eqn(Params,t), 'Color', c2, 'LineWidth', 2)
    xlabel('Time (s)', 'FontSize', 14)
    ylabel({'Avergage height'; '(arb. units)'}, 'FontSize', 14)
    xlim([0 800])
    ylim([0 1.5])
    set(gca, 'XTick', [0 300 600])
    set(gca, 'YTick', [0 0.5 1 1.5])
    set(gca, 'LineWidth', 2, 'FontSize', 14)
    title(['3*\tau = ', num2str(round(tau_3)), ' s'], 'FontSize', 16)
    
    dashline(t0_x_plot, t0_y_plot, 3, 2, 3, 2, 'Color', c3, 'LineWidth', 2)
    dashline(tend_x_plot, tend_y_plot, 3, 2, 3, 2, 'Color', c3, 'LineWidth', 2)

    if save_tau_fits_plots == 1
        fig_name = strcat('SG_filtered_FitAverageHeight_tau_tracknumb_n=', num2str(i));
        saveas(gca, fullfile(outdirectory, fig_name), 'pdf')
        close(gcf);
    end
    
    % save out data and fits as text file
    tau3_x_start_line = ones(length(t), 1) * t0;
    tau3_y_start_line = linspace(0, 10, length(t))';
    tau3_x_end_line   = ones(length(t), 1) * tend;
    tau3_y_end_line   = linspace(0, 10, length(t))';
    
    time_coordinate = t;
    tanh_fit        = Fitting_eqn(Params,t);
    
    pore_formation_height_time_fit = [t, h_vec_sg, tanh_fit, tau3_x_start_line, tau3_y_start_line, tau3_x_end_line, tau3_y_end_line];
    
    data_array_name = strcat('Individual_PoreFormation_PlotArray_n=', num2str(i), '.txt');
    dlmwrite(fullfile(outdirectory, data_array_name), pore_formation_height_time_fit);
    
end

%% Plot scatter - time of pore growth against t_{0}

figure(),
plot(firsttrack_sec_array_final, delta_t_3tau_sg_final, 'x', 'Color', c2, 'LineWidth', 2)
set(gca, 'FontSize', 14, 'LineWidth', 2)
xlabel('t_{0} (s)', 'FontSize', 14)
ylabel('3*\tau (s)', 'FontSize', 14)
title(strcat('3*\tau as function of t_0 ', ' (N=', num2str(length(delta_t_3tau_sg_final)), ')'), 'FontSize', 16)
axis([0 3000 0 500])

%% Create videos of pore growth with mean height pltos and fitting

% again some repetition, but again it works and it's just for making
% figures - so I'm leaving it be!

if save_plotting_videos == 1

    for n = 1:length(mh_track_fitting_final)

        % extract time and height vectors and parameters
        t        = mh_track_fitting_final{n}(:,1);
        h_vec_sg = mh_track_fitting_final{n}(:,2);

        Params        = Params_cell_final{n};
        template_size = template_nm_final(n);
         
        track_seq = img_seq_final{n};

        pause('on')

        Video_Seq_title_savename = strcat(Video_Seq_title_pg, '_', num2str(n));
        Video_Save_FullFileName = fullfile(outdirectory, Video_Seq_title_savename);

        vid = VideoWriter(Video_Save_FullFileName, 'MPEG-4');
        vid.FrameRate = Saved_frames_per_sec;
                
        h = figure();
        open(vid);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.1, 0.3, 0.7]);
        
        for i = 1:length(track_seq)
            hold on
            img    = track_seq{i};
            t_elem = t(i);
            h_elem = h_vec_sg(i);

            subplot(211),
            imagesc(img, 'XData', [0 template_size], 'YData', [0 template_size])
            set(gca, 'YDir', 'norm', 'YTickLabel', [])
            set(gca, 'FontSize', 14)
            colormap('copper')
            caxis(clims)
            pbaspect([1 1 1])
            xlabel('nm', 'FontSize', 14)
            ax2 = subplot(212);
            scatter(t_elem, h_elem, 5, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'LineWidth', 2);
            set(gca, 'FontSize', 14)
            xlim([0 800])
            ax2.XTick = [0, 300, 600];
            ylim([0 1.5])
            ax2.YTick = [0, 0.5, 1, 1.5];
            pbaspect([1 1 1])
            xlabel('Time (s)', 'FontSize', 14)
            ylabel({'Average height'; '(arb. units)'}, 'FontSize', 14)
            
            % to have some time at the end of the video where you can see the result of the fitting
            if i == length(track_seq) 
                for k = 1:25

                    t0   = t0_3tau_sg_final(n);
                    tend = tend_3tau_sg_final(n);

                    t0_x   = [t0 t0];
                    t0_y   =  [0 1.5];
                    tend_x = [tend tend];
                    tend_y =  [0 1.5];
                    
                    t0_x_plot = linspace(t0, t0, 200);
                    t0_y_plot = linspace(0, 1.5, 200);
                    
                    tend_x_plot = linspace(tend, tend, 200);
                    tend_y_plot = linspace(0, 1.5, 200);

                    plot(t, Fitting_eqn(Params,t), 'Color', c2, 'LineWidth', 2)
                    dashline(t0_x_plot, t0_y_plot, 1.8, 1, 1.8, 1, 'Color', c3, 'LineWidth', 2)
                    dashline(tend_x_plot, tend_y_plot, 1.8, 1, 1.8, 1, 'Color', c3, 'LineWidth', 2)

                    Img_frame = getframe(h);
                    writeVideo(vid, Img_frame);

                end

            else
                Img_frame = getframe(h);
                writeVideo(vid, Img_frame);
            end

        end

        hold off
        close(vid);
        close(gcf);

    end

end

%% plot histogram of 3*tau

min_delta_t = 0;
max_delta_t = ceil(max(delta_t_3tau_sg_final)) + binwidth_s;

binedges = min_delta_t:binwidth_s:max_delta_t+30;

% for 3*tau
figure(), 
delta_t_hist_3tau_sg = histogram(delta_t_3tau_sg_final, binedges, 'FaceColor', c3, 'EdgeColor', c3, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 14)
ylabel('Counts', 'FontSize', 14)
title(strcat('Time of individual MAC formation (3*\tau)', ' (N=', num2str(length(delta_t_tau_sg_final)), ')'), 'FontSize', 16)
axis([0 (max_delta_t+30) 0 10])
set(gca, 'FontSize', 13, 'LineWidth', 2)
fig_name = strcat('Hist_n=', num2str(length(mh_track_fitting_final)));
saveas(gca, fullfile(outdirectory, fig_name), 'jpg')

%% Save chosen pores out as images
% if wish to save out the images from one pore

if save_pore_as_image == 1
    
    images = img_seq_final{pore_n};
    
    for i = 1:length(images)
        
        img=images{i};

        figure(),
        imagesc(img, 'XData', [0 template_size], 'YData', [0 template_size])
        set(gca, 'YDir', 'norm', 'YTickLabel', [])
        set(gca, 'FontSize', 14)
        colormap('copper')
        caxis(clims)
        pbaspect([1 1 1])
        xlabel('Distance (nm)', 'FontSize', 14)

        fig_name = strcat('Pore_numb=', num2str(n), '_img_numb=', num2str(i));
        saveas(gca, fullfile(outdirectoryimgseq, fig_name), 'jpg')

        close(gcf)

    end
end

display(['Mean time of pore formation = ', num2str(mean(delta_t_3tau_sg_final)), ' s, from N = ', (num2str(length(img_seq_final))), ' events'])

