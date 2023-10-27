function array = unjumble(cells, thresh)
    % Aim is to take the cell array that is produced when solving and
    % convert it into something usable. Problem is that one column isnt
    % necessarily one mode, so some kind fo algorithm has to be used to
    % figure out where the numbers go in each row. Thresh is the max jump
    % in the phase velocity permitted in a single mode.

    % figure out how bug the array needs to be
    Ncols = 0;
    for ii = 1:length(cells)
        if length(cells{ii})>Ncols
            Ncols = length(cells{ii});
        end
    end
    array = nan(length(cells), Ncols);
    
    % place first row
    array(1,1:length(cells{1})) = cells{1};

    % loop through rows
    for ii = 2:length(cells)
        to_place = sort(cells{ii}); % values to place in the array, increasing order
        placed = false(size(to_place)); % bool to keep track of which have placed
        compVals = array(ii-1,:); % values to compare to, start with prev row
        for jj = 1:length(compVals)
            if isnan(compVals(jj))
                % look for a non-nan value above it to replace with in case
                % there was a gap
                column = array(:,jj);
                if any(~isnan(column)) % there is one
                    compVals(jj) = column( find(~isnan(column), 1, 'last') );
                end
            end
        end

        % loop through values to place
        for jj = 1:length(to_place)
            [delta, col] = min(abs(to_place(jj)-compVals)); % difference between value and prev row
            if delta<thresh
                array(ii,col) = to_place(jj);
                placed(jj) = 1;
            end
        end
        % remove placed values
        to_place = to_place(~placed);
        if ~isempty(to_place)
            % find a completely empty column
            start = findNewCol(array);
            % place remaining in leftover slots
            array(ii,start:start+length(to_place)-1) = to_place(:);
        end
    end

end
function col = findNewCol(array)
    % Need to find the leftmost column which contains only nans for when we
    % need a new column
    bool = false(size(array(1,:))); % one el for each col
    for ii = 1:length(bool)
        bool(ii) = all(isnan(array(:,ii))); % only 1 if all rows are nan in that col
    end
    col = find(bool, 1, 'first'); % return first empty col
end