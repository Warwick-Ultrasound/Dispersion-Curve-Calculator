% Have a nan problem - if there are nans before a jump it's hard to detect.

clear;
clc;
close all;

filename = "All_modes_12_05_23.csv";

% threshold for jump detection
thresh = 300; % in m/s
smallestSegment = 20; % smallest section to try to home

% read in data
data = readmatrix(filename);
freq = data(:,1);
modes = data(:,2:end);
clear data;

% look for starts and stops
starts = cell(size(modes, 2), 1);
stops = starts;
for mm = 1:size(modes, 2)

    mode = modes(:,mm);
    modeStart = find(~isnan(mode),1); % first non-nan entry

    starts{mm} = [];
    stops{mm} = [];
    
    for ii = modeStart:length(mode)-1
        
        % if mode(ii) is nan, find values either side
        if isnan(mode(ii))
            leftI = modeStart - 1 + find(~isnan(mode(modeStart:ii)), 1, 'last'); % last non-nan before ii
            rightI = ii + find(~isnan(mode(ii:end)), 1, 'first'); % first non-nan after ii
            left = mode(leftI);
            right = mode(rightI);
        else
            leftI = ii;
            rightI = ii+1;
            left = mode(leftI);
            right = mode(rightI);
        end
        
        if left - right > thresh
            stops{mm}(end+1) = rightI-1;
        elseif right - left > thresh
            starts{mm}(end+1) = leftI+1;
        end
    end
end

% remove duplicates cause when we are in a region of nans
for mm = 1:size(modes, 2)
    starts{mm} = unique(starts{mm});
    stops{mm} = unique(stops{mm});
end

% now filter out any tiny bits
delete_list = cell(size(modes, 2), 1);
for mm = 1:size(modes, 2)
   delete_list{mm} = [];
   if ~isempty(starts{mm})
       for ii = 1:length(starts{mm})
           if stops{mm}(ii)-starts{mm}(ii) < smallestSegment
               delete_list{mm}(end+1) = ii;
               % remove segment from modes
               modes(starts{mm}(ii):stops{mm}(ii), mm) = 0/0;
           end
       end
   end
end

for mm = 1:size(modes, 2)
    if ~isempty(delete_list{mm})
        starts{mm}(delete_list{mm}) = [];
        stops{mm}(delete_list{mm}) = [];
    end
end

% now pull out all the lost segments ready to rehome them
lost_segments = cell(size(modes, 2), 1);
for mm = 1:size(modes, 2)
    for ii = 1:length(starts{mm})
        lost_segments{mm, ii} = modes(starts{mm}(ii):stops{mm}(ii),mm);
        lost_segments{mm,ii}(isnan(lost_segments{mm,ii})) = []; % remove any nans 
        modes(starts{mm}(ii):stops{mm}(ii), mm) = 0/0;
    end
end

% check that they seem reasonable
figure;
tiledlayout(4,4);
for ii = 1:size(modes, 2)
    nexttile;
    plot(modes(:,ii));
    hold on;
    if ~isempty(starts{ii})
        for jj = 1:length(starts{ii})
            xregion(starts{ii}(jj), stops{ii}(jj));
        end
    end
    hold off;
    title(ii)
end

% find all the starts and stops of the gaps
gapStarts = cell(size(modes, 2), 1);
gapStops = gapStarts;
for mm = 1:size(modes, 2)
    gapStarts{mm} = [];
    gapStops{mm} = [];
    modeStart = find(~isnan(modes(:,mm)), 1); % start of trace

    for ii = modeStart:length(modes(:,mm))-1
        if ~isnan(modes(ii,mm)) && isnan(modes(ii+1,mm))
            gapStarts{mm}(end+1) = ii+1;
        elseif isnan(modes(ii,mm)) && ~isnan(modes(ii+1,mm))
            gapStops{mm}(end+1) = ii;
        end
    end
end

% delete tiny gaps
del_list = cell(size(modes, 2), 1);
for mm = 1:size(modes, 2)
    del_list{mm} = [];
    for ii = 1:length(gapStarts{mm})
        if gapStops{mm}(ii) - gapStarts{mm}(ii) < smallestSegment
            del_list{mm}(end+1) = ii;
        end
    end
end

for mm = 1:size(modes, 2)
    if isempty(del_list{mm})
        continue;
    end
    gapStarts{mm}(del_list{mm}) = [];
    gapStops{mm}(del_list{mm}) = [];
end

% now loop through the lost segments and find which mode gaps they might
% fit into best
for mseg = 1:size(modes, 2) % loop through modes the segments came from
    for iseg = 1:size(lost_segments, 2) % segment number
        
        % here is where we find a home for that segment
        seg = lost_segments{mseg, iseg};
        
        if isempty(seg) % dont do anything if empty, no gaps
            continue;
        end

        % loop through the gaps
        for mgap = 1:size(modes, 2)
            searchMode = modes(:,mgap);
            for igap = 1:length(gapStarts{mgap})
                ii = gapStarts{mgap}(igap); % start of gap
                jj = gapStops{mgap}(igap); % end of gap

                if abs(mean(seg) - mean([searchMode(ii-1), searchMode(jj+1)])) < 20
                    disp("Found a good match!");
                    disp("from mode "+string(mseg));
                    disp("Going into mode "+string(mgap));
                    disp("length of gap is "+string(jj-ii));
                    disp("length of segment is "+string(length(seg)));
                    disp(' ');

                    modes(ii:ii+length(seg)-1, mgap) = seg;
                end

            end
        end

    end
end

% check that they seem reasonable
figure;
tiledlayout(4,4);
for ii = 1:size(modes, 2)
    nexttile;
    plot(modes(:,ii));
    hold on;
    if ~isempty(starts{ii})
        for jj = 1:length(starts{ii})
            xregion(starts{ii}(jj), stops{ii}(jj));
        end
    end
    hold off;
    title(ii)
end

% now we just need to fill in the nans with spline
for mm = 1:size(modes, 2)
    startMode = find(~isnan(modes(:,mm)), 1);
    mode = modes(startMode:end, mm);
    freqCrop = freq(startMode:end);

    mask = isnan(mode);
    freqCropGap = freqCrop;
    freqCropGap(mask) = [];
    mode(mask) = [];
    newMode = spline(freqCropGap, mode, freqCrop);

    if sum(isnan(newMode)) > 0
        disp(mm);
    end

    modes(startMode:end, mm) = newMode;
end

% check that they seem reasonable
figure;
tiledlayout(4,4);
for ii = 1:size(modes, 2)
    nexttile;
    plot(modes(:,ii));
    hold on;
    if ~isempty(starts{ii})
        for jj = 1:length(starts{ii})
            xregion(starts{ii}(jj), stops{ii}(jj));
        end
    end
    hold off;
    title(ii)
end

% save the data
savedata = [freq, modes];
outfile = split(filename, '.');
outfile = outfile(1) + "GAP_FILLED.csv";
writematrix(savedata, outfile);