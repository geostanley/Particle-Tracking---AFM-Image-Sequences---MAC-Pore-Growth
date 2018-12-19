%%%=== XYZarray_indexed_by_percentage_height ===%%%

% This function takes in an Nx3 array (i.e., XYZ coordinates), takes a
% histogram of the data, and only keeps the data which is either above
% (greater_than == 1), or below (greater_than == 0) the Planefit_mask,
% which is a number between 0 and 1, and is a percentage. 

% User must also input the BinWidth in nm. The smaller the better. Value
% should probably be <0.5 nm.

% I.e., greater_than = 1 & Plane_fit_mask = 0.4, means keep the top 60% of
% the height data.

function [XYZ_array_for_plane_fit, XYZ_array_not_for_plane_fit] = XYZarray_indexed_by_percentage_height(XYZ_array, BinWidth_nm, Plane_fit_mask, greater_than)

    X_array = XYZ_array(:,1);
    Y_array = XYZ_array(:,2);
    Z_array = XYZ_array(:,3);

    % Find the top x% of height data for plane fitting.
    % Take a histogram of data, find the desired cut-off point, and only take
    % that data forward for finding the 1st order plane fit - just like
    % Nanoscope Analysis.

    % take histogram of height data, binning into 1nm bins
    Z_histedges = [floor(min(Z_array)):BinWidth_nm:ceil(max(Z_array))];
    Z_hist_counts = histcounts(Z_array, Z_histedges);
    Z_hist_counts_sum = sum(Z_hist_counts(:));
    
    % create an array of bin centres in nm
    Bin_centres = zeros(1, length(Z_histedges)-1);
    for i = 1:length(Bin_centres)
        Bin_centres(i) = Z_histedges(i) + (BinWidth_nm/2);
    end

    % create an array with aggregate counts based on bins
    Z_aggregate_hist_counts = zeros(size(Z_hist_counts));
    Z_aggregate_hist_counts(1) = Z_hist_counts(1);
    for i=1:length(Z_hist_counts)-1
        Z_aggregate_hist_counts(i+1) = Z_aggregate_hist_counts(i) + Z_hist_counts(i+1);
    end

    % represent aggregate histcounts as running percent, then find the 
    % cut-off point as an index and in nm
    Z_aggregate_hist_counts_percentage = Z_aggregate_hist_counts ./ Z_hist_counts_sum;
    Z_aggregate_hist_counts_abs = abs(Z_aggregate_hist_counts_percentage - Plane_fit_mask);
    Z_cut_off_idx = find(Z_aggregate_hist_counts_abs == min(Z_aggregate_hist_counts_abs),1);  
    Z_cut_off_nm  = Bin_centres(Z_cut_off_idx);

    % Pull out the data points at, or above the cut off point (nm)
    if greater_than == 1
        Z_array_for_plane_fit_idx = find(Z_array >= Z_cut_off_nm);
        Z_array_not_for_plane_fit_idx = find(Z_array < Z_cut_off_nm);
    else
        Z_array_for_plane_fit_idx = find(Z_array < Z_cut_off_nm);
        Z_array_not_for_plane_fit_idx = find(Z_array >= Z_cut_off_nm);
    end

    % save desired data
    X_array_for_plane_fit = X_array(Z_array_for_plane_fit_idx);
    Y_array_for_plane_fit = Y_array(Z_array_for_plane_fit_idx);
    Z_array_for_plane_fit = Z_array(Z_array_for_plane_fit_idx);
    XYZ_array_for_plane_fit = [X_array_for_plane_fit Y_array_for_plane_fit Z_array_for_plane_fit];
    
    % save the rest of the data for plotting if wanted
    X_array_not_for_plane_fit = X_array(Z_array_not_for_plane_fit_idx);
    Y_array_not_for_plane_fit = Y_array(Z_array_not_for_plane_fit_idx);
    Z_array_not_for_plane_fit = Z_array(Z_array_not_for_plane_fit_idx);
    XYZ_array_not_for_plane_fit = [X_array_not_for_plane_fit Y_array_not_for_plane_fit Z_array_not_for_plane_fit];
    
end