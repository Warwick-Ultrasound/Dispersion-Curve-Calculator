

% ============================================================
% --------------------- Code description ---------------------
% Dispersion curves according to Analytical approach in Copper
% plate samples of which the geometrical characteristics and
% the Longitudinal and Shear ultrasonic velocities are
% determined experimentally . This code is executed to
% calculate the phase velocity ,c_p , over arbitrary defined
% frequency and phase velocity range with an increasing steps,
% respectively . Successively, the remaining needed values of
% circular wavenumber and group velocit are determined. These
% are used to plot the respective dispersion curve diagrams.
% ============================================================
clear;
clc;
close all;

saveFile = 'OutputFileNameHere.csv';

tic; % To count the time needed to execute the code .

T = 15; % temperature ( delete and replace c_long and c_shear if required )

c_long = c_cu_long(T); % Longitudinal velocity in material (m/s)
c_shear = c_cu_shear(T); % Shear velocity in material (m/s)
thickness = 0.7E-3; % Full plate thickness /m

thick = thickness/2; % half plate thickness of the plate, as this is what the equations use (h)

freq = (0:1000:5e6); % Frequency range (Hz)
c_phase = (0.1:10:8000.0); % Phase velocity range (m/s)

% ============================================================


% We loop the process over the two ranges. Here we define the symmetric (S)
% and anti-symmetric (A) modes. The inner loop calculates these modes for
% increments of c_phase (the range of phase velocities). The outer loop
% runs the inner loop for incrementing values of freq (the frequency
% range). The inner loop calculates in three regions: c_p < c_shear, c_p
% between c_shear and c_long, and c_p > c_long. Equations are given in
% Chapter 8 of J Rose 'Ultrasonic Waves in Solid Media'

circ_freq = freq*2*pi;
k = zeros(length(c_phase), length(freq));
p = k;
q = k;
symr = k;
anti_symr = k;

for i = 1:length(freq)
    disp(string(i/length(freq)*100)+"% done");
    w = circ_freq(i);

    parfor j = 1:length(c_phase)

        k(j,i)= w/c_phase(j) ;   %calculate wave number

        % First Region
        if c_phase(j)<c_shear

            
            p(j,i)= 1i*w* sqrt(((c_long^2) - (c_phase(j).^2))/((c_phase(j).^2)*(c_long^2)));
            q(j,i)= 1i*w*sqrt(((c_shear^2) - (c_phase(j).^2))/((c_phase(j).^2)*(c_shear^2)));
            symr(j,i)=((tan(q(j,i).*thick))/q(j,i))+(4*(k(j,i).^2)*p(j,i).*tan(p(j,i).*thick))/(((q(j,i).^2)-(k(j,i).^2))^2);
            anti_symr(j,i)=(q(j,i).*tan(q(j,i).*thick))+(((q(j,i).^2-k(j,i).^2)^2)*tan(p(j,i).*thick))/(k(j,i).^2*p(j,i).*4);

            % Second Region
        elseif c_phase(j)<c_long && c_phase(j)>c_shear

            p(j,i)= 1i*w* sqrt(((c_long^2) - (c_phase(j).^2))/((c_phase(j).^2)*(c_long^2)));
            q(j,i)= sqrt(((w/c_shear)^2)-((w/c_phase(j))^2));
            symr(j,i)=((tan(q(j,i).*thick))/q(j,i))+(4*(k(j,i).^2)*p(j,i).*tan(p(j,i).*thick))/(((q(j,i).^2)-(k(j,i).^2))^2);
            anti_symr(j,i)=(q(j,i).*tan(q(j,i).*thick))+((((q(j,i).^2)-(k(j,i).^2))^2)*tan(p(j,i).*thick))/((k(j,i).^2)*p(j,i).*4);

            % Third Region
        else

            p(j,i)=sqrt(((w/c_long)^2)-((w/c_phase(j))^2));
            q(j,i)=sqrt(((w/c_shear)^2)-((w/c_phase(j))^2));
            symr(j,i)=((tan(q(j,i).*thick))/q(j,i))+(4*(k(j,i).^2)*p(j,i).*tan(p(j,i).*thick))/(((q(j,i).^2)-(k(j,i).^2))^2);
            anti_symr(j,i)=(q(j,i).*tan(q(j,i).*thick))+((((q(j,i).^2)-(k(j,i).^2))^2)*tan(p(j,i).*thick))/((k(j,i).^2)*p(j,i).*4);

        end
    end
end
% ============================================================



%Phase velocity calculation
%Locates points where the symmetric and anti-symmetric modes cross x axis -
%these are the values of phase velocity (C_p)

%preallocate arrays
symr_all = zeros(length(c_phase),length(freq));
antisymr_all = zeros(length(c_phase),length(freq));

for ii = 1:length(freq)

    A = symr (:,ii);
    B = anti_symr(:,ii);

    for cc = 1:(length(c_phase)-2)

        if A(cc,:) < 0 && A(cc+1,:) > 0 && A(cc+2,:)>A(cc+1,:) && A(cc-1,:)< A(cc,:)
            symr_all(cc,ii) = cc;
        elseif A(cc,:) > 0 && A(cc+1,:) < 0 && A(cc+2,:) < A(cc+1,:) && A(cc-1,:)> A(cc,:)
            symr_all(cc,ii) = cc;
        else
            symr_all(cc,ii) = 0;
        end

        if B(cc,:) < 0 && B(cc+1,:) >0 && B(cc+2,:)>B(cc+1,:) && B(cc-1,:)< B(cc,:)
            antisymr_all(cc,ii) = cc;
        elseif B(cc,:) > 0 && B(cc+1,:) <0 && B(cc+2,:)< B(cc+1,:) && B(cc-1,:)> B(cc,:)
            antisymr_all(cc,ii) = cc;
        else
            antisymr_all(cc,ii) = 0;
        end

    end
end

zz = cell(1,length(freq));
yy = cell(size(zz));

for ii = 1:length(freq)
    zz{ii} = find(symr_all(:,ii) ~= 0);
    yy{ii} = find(antisymr_all(:,ii) ~= 0);
end


x_symr = NaN( length(zz{ii}),length(freq));
x_antisymr = NaN(length(yy{ii}),length(freq));


for jj = 1:length(freq)
    for ii = 1:length(zz{jj})
        x_symr(ii,jj)= zz{jj}(ii);
    end
    for ii = 1:length(yy{jj})
        x_antisymr(ii,jj)= yy{jj}(ii);
    end
end

% display execution time
time = toc;
disp(string(time)+" seconds taken to calculate");

M = freq.*(thick*2);

x_symr_1 = x_symr;
x_symr_1(isnan(x_symr_1))=1;

cp_symr_1 = zeros(size(x_symr_1,1), length(freq));

for ss = 1:size(x_symr_1,1)
    for cc = 1:length(freq)
        cp_symr_1(ss,cc) = c_phase(x_symr_1(ss,cc));
    end
end
cp_symr_1(cp_symr_1 == c_phase(1)) = NaN;

x_antisymr_1 = x_antisymr;
x_antisymr_1(isnan(x_antisymr_1))=1;
cp_antisymr_1 = zeros(size(x_antisymr_1,1), length(freq));
for sa = 1:size(x_antisymr_1,1)
    for cc = 1:length(freq)
        cp_antisymr_1(sa,cc) = c_phase(x_antisymr_1(sa,cc));
    end
end
cp_antisymr_1(cp_antisymr_1 == c_phase(1)) = NaN;

figure()
hold on
plot(freq,cp_symr_1,'r')
plot(freq,cp_antisymr_1,'b')
legend('Symmetric','Antisymmetric','Location','northwest');
ylabel('Phase velocity')
xlabel('Frequency')
title('1 mm Al plate Lamb Dispersion curve')

fd = freq./1e6 * thick/1E-3; % freq thickness in MHz.mm

fd_save = fd';
A_save = cp_antisymr_1';
S_save = cp_symr_1';

outdata = [fd_save, S_save, A_save];

disp("SPLIT = "+string(1+size(S_save, 2))); % split is the column which separates sym and asym modes, need to input to calc_Cg

writematrix(outdata, saveFile);