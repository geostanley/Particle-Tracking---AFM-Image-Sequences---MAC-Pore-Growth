function appearancetable = createappearancetable(trackingtable,frames)
% Based on an input trackingtable with x coord, y coord, frame no, tracking
% ID in its 1st, 2nd, 3rd and 4th column, this function selects the first
% "postframes" track entries counting from the first frame in which a track
% appears, and adds entries (with same coordinates) for the last
% "preframes" prior to this appearance.
% The frames input parameter should be of the form [preframes postframes].
% It assumes that the tracking IDs are positive integer numbers.
% Output is in the same format at trackingtable.

preframes  = frames(1);
postframes = frames(2);

% Sort input array first according to 4th column (track ID) and next to 3th
% column (frame number)
sortedtable = sortrows(trackingtable,[4 3]);

% Fill in gaps, using values for next frame in case of missed frames
difftable = diff(sortedtable(:,3)); % Calculate differences of frame numbers
indices   = find(difftable>1)+1; % Check for frame gaps

while ~isempty(indices)
    
    inserttable      = sortedtable(indices,:); % Read rows after "frame gaps"
    inserttable(:,3) = inserttable(:,3)-1; % Change frame number by -1
    sortedtable      = insertrows(sortedtable,inserttable,indices-1); % Insert entries (identifical to subsequent entry)
    difftable        = diff(sortedtable(:,3)); % Recalculate differences of frame numbers
    indices          = find(difftable>1)+1; % Re check for single-frame gaps
    
end

% Create output table
appearancetable=[];

% Run over tracks
for i=1:max(sortedtable(:,4))
    
    indices = find(sortedtable(:,4)==i); % Find entries for track i
    indices = indices(1:min(length(indices),postframes)); % Only consider first postframes entries
    
    if (~isempty(indices)) %& (sortedtable(indices(1),3)>1)
        firstentry = sortedtable(indices(1),:); % First entry
        firstframe = firstentry(1,3);
        
        for j = max(1,firstframe-preframes):(firstframe-1)
            priorentry      = firstentry;
            priorentry(1,3) = j; % Set frame numbers before appearance
            % Next add entries before appearance
            appearancetable = cat(1,appearancetable,priorentry);
        end
        % And then add first postframe entries after appearance
        appearancetable=cat(1,appearancetable,sortedtable(indices,:));
        
    end
end

