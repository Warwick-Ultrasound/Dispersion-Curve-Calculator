function [Sroots, ASroots] = solver(w, h, cl, ct, cp_range)
    % A function which finds roots of symm and asymm within a range in cp
    % at a particular value of angular frequency w. Inputs are angular
    % frequency, half thickness of plate, longitudinal speed of sound,
    % shear speed of sound, and a 2-vector with the range of cp values to
    % check.

    opt = optimset('fzero'); % solver options
    opt.TolX = 1E-30; % lower tolerance than default
    opt.FunValCheck = 'on'; % check for NaNs, infs etc when solving
    opt.Display = 'off'; % supress message if no root found

    % generate an array of cp values to test initially
    cp_list = linspace(cp_range(1), cp_range(2), 10000);

    % solve for symmetric roots first
    Sval = symm(cp_list, w, h, cl, ct); % y values of symm fucntion, find roots in this
    % find approx roots first
    Sroots = [];
    for ii = 2:length(Sval)-2
        upcross = Sval(ii) < 0 && Sval(ii+1) > 0 && Sval(ii+2)>Sval(ii+1) && Sval(ii-1)<Sval(ii);
        downcross = Sval(ii) > 0 && Sval(ii+1) < 0 && Sval(ii+2)<Sval(ii+1) && Sval(ii-1)>Sval(ii);
        if upcross || downcross % zero corssing
            Sroots(end+1) = cp_list(ii);
        end
    end
    
    for ii = 1:length(Sroots)
         [root, ~, flag] = fzero(@symm, Sroots(ii), opt, w, h, cl, ct); % recalculate to better tolerance
         if flag == 1
            Sroots(ii) = root;
         end
    end

    % now asymmetric
    ASval = asymm(cp_list, w, h, cl, ct); % y values of symm fucntion, find roots in this
    % find approx roots first
    ASroots = [];
    for ii = 2:length(ASval)-2
        upcross = ASval(ii) < 0 && ASval(ii+1) > 0 && ASval(ii+2)>ASval(ii+1) && ASval(ii-1)<ASval(ii);
        downcross = ASval(ii) > 0 && ASval(ii+1) < 0 && ASval(ii+2)<ASval(ii+1) && ASval(ii-1)>ASval(ii);
        if upcross || downcross % zero corssing
            ASroots(end+1) = cp_list(ii);
        end
    end
    for ii = 1:length(ASroots)
         [root, ~, flag] = fzero(@asymm, ASroots(ii), opt, w, h, cl, ct); % recalculate to better tolerance
         if flag == 1
            ASroots(ii) = root;
         end
    end
end