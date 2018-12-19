function trackingtable = createtracktable(coordinates,tracks)
% Put track results from simpletracker.m put into matrix with 
% col1 x; col2 y; col3 framenumber; col4 trackID ...
% Inputs are:
% <coordinates>: A cell array with each cell corresponding to an image
% frame, anc containing a list of coordinates of features detected in that
% particular frame
% <tracks>: A cell array with each cell corresponding to a track,
% and containing an array of the same length as coordinates, with the ith
% element giving the index of the track point in coordinates{i}

trackingtable=[];
for i=1:length(tracks)
    keeptrack=tracks{i}; % Read track i
    coordinateindices=find(~isnan(keeptrack)); % Find not NaN entries
    for j=1:length(coordinateindices)
        keepcoordinates=coordinates{coordinateindices(j)}; % Read for frame
        if length(keepcoordinates)>0
            rowelements=cat(2,keepcoordinates(keeptrack(coordinateindices(j)),:),coordinateindices(j),i);
            trackingtable=cat(1,trackingtable,rowelements);
        end;
    end;
end;

end

