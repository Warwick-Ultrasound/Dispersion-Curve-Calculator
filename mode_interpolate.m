function modes = mode_interpolate(freq, modes)
    %  The task for this function is to find the gaps in the modes where
    %  there are nans and fill them in using spline interpolation

    for mm = 1:size(modes, 2)   % cycle through modes
        
        mode = modes(:,mm);

        start = find(~isnan(mode), 1, 'first'); % find first non-nan value

        freq_crop = freq(start:end); % from start -> end
        freq_crop2 = freq_crop; % will have gaps removed
        mode_crop = mode(start:end);
        if length(mode_crop)>201
            mode_crop = sgolayfilt(mode_crop, 2, 201);
        end
        nanBool = isnan(mode_crop); % indicates where gaps are
        freq_crop2(nanBool) = []; % remove freq values where mode is nan
        mode_crop(nanBool) = []; % remove mode values where mode is nan
        fixed_mode = spline(freq_crop2, mode_crop, freq_crop); % fill in the gaps

        modes(start:end,mm) = fixed_mode; % put it back in
        
    end
end