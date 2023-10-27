function group_modes = groupVel(freq, full_thick, modes)
    % Takes the phase velocities and converts them to group velocities
    % using the equation in Rose Ultrasonic Waves in Solid Media.

    fd = full_thick*freq;
    fd_crop = fd(1:end-1);
    d_fd = fd(2)-fd(1);

    group_modes = zeros([size(modes, 1)-1, size(modes, 2)]); % cut a row off due to diff

    for ii = 1:size(modes, 2)
        mode = modes(:,ii);
        mode = sgolayfilt(mode, 2, 201);
        diff_mode = diff(mode)/d_fd;
        diff_mode = sgolayfilt(diff_mode, 2, 201);
        mode = mode(1:end-1); % same length as diff_mode
        group_modes(:,ii) = mode.^2.*((mode - fd_crop.*(diff_mode)).^-1);
    end

end