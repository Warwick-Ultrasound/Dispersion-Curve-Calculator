% reads in an array which contains the frequency in column one, followed by
% all of the modes phase velocities in columns 2:end. Converts to a group
% velocity curve.

clear;
clc;
close
cp_data = readmatrix("iriginal_file_name_GAP_FILLED.csv");
freq = cp_data(:,1);
cp_data = cp_data(:,2:end);
circ_freq = 2*pi*freq;
thickness = 0.7E-3;
fd = (freq/1E6)*(thickness/1E-3); % in MHz mm

figure;
plot(fd,cp_data)
title('Dispersion curves, phase velocity (Cp)');
xlabel('Frequency-thickness, MHz-mm');
ylabel('Phase Velocity, m/s)');

% if you look closely the input has steps in it, which causes derivatives
% problems => smooth

% smoothing parameters
order = 3;
winLen = 41;
for ii = 1:size(cp_data, 2)
    cp_data(:,ii) = sgolayfilt(cp_data(:,ii), order, winLen);
end

%------------CALCULATE GROUP VELOCITIES----------------------------%

d_fd = fd(2) - fd(1);
cg_modes = zeros(size(cp_data, 1)-1, size(cp_data, 2)); % shorter due to diff
fd_crop = fd(1:end-1);
figure;
hold on;
for ii = 1:size(cp_data, 2)
    mode = cp_data(:,ii);
    diff_mode = diff(mode)/d_fd;
    mode = mode(1:end-1); % so its the same length as diff

    diff_mode = sgolayfilt(diff_mode, 3, 401);

    cg_mode = mode.^2.*((mode - fd_crop.*(diff_mode)).^-1);

    cg_mode = sgolayfilt(cg_mode, 1, 351);

    plot(abs(fftshift(fft(cg_mode))));
    cg_modes(:,ii) = cg_mode(:);
end
hold off;

split = 7;
sym = cg_modes(:,1:split);
asym = cg_modes(:,split+1:end);

figure;
plot(fd_crop,sym, 'b-', fd_crop, asym, 'r-');
title('Dispersion curves, Group velocity');
xlabel('Frequency-thickness, MHz-mm');
ylabel('Group Velocity, m/s)');
xlim([0 8]);
