%{
The main script that uses the solver to create dispersion curves.
%}

clear;
clc;
close all;

% INPUTS

cp_min = 100; % minimum phase velocity to calculate, make larger if getting negative cp
cp_max = 10E3; % max phase velocity to calculate
f_min = 100E3; % minimum frequency
f_max = 10E6; % maximum frequency
f_step = 1E3; % step size between calculation frequencies

c_L = 4642; % Longitudinal sound speed.
c_S = 2389; % Shear sound speed. 
d = 1E-3; % Full plate thickness (1mm if you want frequency thickness product on x)

material = "Copper"; % material for generating titles and filenames

% CALCULATION BEGIN
cp_range = [cp_min, cp_max];
f = (f_min:f_step:f_max).';
w = 2*pi*f;
h = d/2; % half plate thickness

% begin calculation loop for c_p
Scell = cell(length(w), 1);
AScell = cell(length(w), 1);
parfor ii = 1:length(w)
    [Scell{ii}, AScell{ii}] = solver(w(ii), h, c_L, c_S, cp_range);
end

% re-jig results into a proper array to plot
S = unjumble(Scell, 600); % number is biggest jump in c_p allowed from one value to the next in a particular mode
AS = unjumble(AScell, 600);

% plot results
figure;
for ii = 1:size(AS, 2)
    plot(f/1E6, AS(:,ii), 'DisplayName', "A"+string(ii-1));
    if ii == 1
        hold on;
    end
end
for ii = 1:size(S, 2)
    plot(f/1E6, S(:,ii), 'DisplayName', "S"+string(ii-1));
end
hold off;
ylim(cp_range);
legend;
title("Check Mode Categorisation", "Each mode should be its own colour, adjust line 37,38");
xlabel('Frequency /MHz');
ylabel('Phase Velocity /ms^{-1}');

% Fill in any small gaps
S = mode_interpolate(f, S);
AS = mode_interpolate(f, AS);

% plot results
figure;
for ii = 1:size(AS, 2)
    plot(f/1E6, AS(:,ii), '-','DisplayName', "A"+string(ii-1));
    if ii == 1
        hold on;
    end
end
for ii = 1:size(S, 2)
    plot(f/1E6, S(:,ii), '--','DisplayName', "S"+string(ii-1));
end
hold off;
ylim(cp_range);
legend;
title(material+" Lamb Wave Dispersion Curves");
xlabel('Frequency /MHz');
ylabel('Phase Velocity /ms^{-1}');

% Convert phase velocities to group velocities
Sgroup = groupVel(f, d, S);
ASgroup = groupVel(f, d, AS);

% plot group velcoities
figure;
for ii = 1:size(ASgroup, 2)
    plot(f(1:end-1)/1E6, ASgroup(:,ii), '-','DisplayName', "A"+string(ii-1));
    if ii == 1
        hold on;
    end
end
for ii = 1:size(Sgroup, 2)
    plot(f(1:end-1)/1E6, Sgroup(:,ii), '--','DisplayName', "S"+string(ii-1));
end
hold off;
legend;
title(material+" Lamb Wave Dispersion Curves");
xlabel('Frequency /MHz');
ylabel('Group Velocity /ms^{-1}');

% ---- save data to ouput file ----

% collate data
phaseData = array2table([f, AS, S]);
groupData = array2table([f(1:end-1), ASgroup, Sgroup]);

% create header names for the columns
NS = size(S, 2); % number symmetric modes
NAS = size(AS, 2); % Number AS modes
cols = cell(1, size(phaseData, 2));
cols{1} = 'Frequency';
for ii = 1:NAS
    cols{ii+1} = char('A'+string(ii-1));
end
for ii = 1:NS
    cols{1+NAS+ii} = char('S'+string(ii-1));
end
phaseData.Properties.VariableNames = cols;
groupData.Properties.VariableNames = cols;

% generate filenames
baseName = material+"_lamb_disp_"+string(d/1E-3)+"mm_";
phaseName = baseName + "phase_vel.csv";
groupName = baseName + "group_vel.csv";

% write to file
writetable(phaseData, phaseName);
writetable(groupData, groupName);
disp("Data written to files.")