% This function finds peaks of cross-correlation of the image
% with a template, and returns the coordinates of cross-correlation
% peaks above the threshold variable.
% It also plots the results in figure(1)

% This function uses FastPeakFind.m (downloaded from MathWorks file exchange)

function coordinates = findfeatures_spmFiles(image,template,threshold)

    % Cross-correlate with template
    C = normxcorr2(template, image);

    % Remove padding added by cross-correlation
    C(1:size(template,1)/2,:)=[]; C(:,1:size(template,2)/2)=[];

    % Set cross-correlation values below "threshold" value to zero
    C(find(C<threshold)) = 0;

    % Use FastPeakFind.m routine (downloaded from MathWorks file exchange)
    p=FastPeakFind(C);

    % Rearrange to have list of (x,y) coordinates for maxima in cross-correlation
    if length(p)>0 p=[p(1:2:end),p(2:2:end)]; end;
    coordinates=p; % Store coordinates

end
