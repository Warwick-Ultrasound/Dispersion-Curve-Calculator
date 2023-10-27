function val = asymm(cp, w, h, cl, ct)
    % A function which evaluates the LHS of eqn 8.38 from Rose.

    k = w./cp;
    kL = w/cl;
    kT = w/ct;
    
    val = zeros(size(k));

    for ii = 1:length(k)
        if k(ii)>kT % REGION I

            p = 1i*w* sqrt(((cl^2) - (cp(ii)^2))/((cp(ii)^2)*(cl^2)));
            q = 1i*w*sqrt(((ct^2) - (cp(ii).^2))/((cp(ii)^2)*(ct^2)));
            val(ii) = q*tan(q*h)+(((q^2-k(ii)^2)^2)*tan(p*h))/(k(ii)^2*p*4);

        elseif k(ii) < kT && k(ii) > kL % REGION II

            p = 1i*w* sqrt(((cl^2) - (cp(ii).^2))/((cp(ii).^2)*(cl^2)));
            q = sqrt(((w/ct)^2)-((w/cp(ii))^2));
            val(ii) = (q*tan(q*h))+((((q^2)-(k(ii)^2))^2)*tan(p*h))/((k(ii)^2)*p*4);
            
        else % REGION III
            
            p = sqrt(((w/cl)^2)-((w/cp(ii))^2));
            q = sqrt(((w/ct)^2)-((w/cp(ii))^2));
            val(ii) = q*tan(q*h) + ((((q^2)-(k(ii)^2))^2)*tan(p*h))/((k(ii)^2)*p*4);

        end
    end
end