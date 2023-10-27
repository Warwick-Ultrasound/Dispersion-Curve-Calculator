function val = symm(cp, w, h, cl, ct)
    % A function which evaluates the LHS of eqn 8.38 from Rose.

    k = w./cp;
    kL = w/cl;
    kT = w/ct;
    
    val = zeros(size(k));

    for ii = 1:length(k)
        if k(ii)>kT % REGION I

            p = 1i*w* sqrt(((cl^2) - (cp(ii)^2))/((cp(ii)^2)*(cl^2)));
            q = 1i*w*sqrt(((ct^2) - (cp(ii).^2))/((cp(ii)^2)*(ct^2)));
            val(ii) = ((tan(q*h))/q)+(4*(k(ii)^2)*p.*tan(p.*h))/(((q^2)-(k(ii)^2))^2);

        elseif k(ii) < kT && k(ii) > kL % REGION II

            p = 1i*w* sqrt(((cl^2) - (cp(ii).^2))/((cp(ii).^2)*(cl^2)));
            q = sqrt(((w/ct)^2)-((w/cp(ii))^2));
            val(ii) = ((tan(q*h))/q)+(4*(k(ii)^2)*p.*tan(p*h))/(((q^2)-(k(ii)^2))^2);
            
        else % REGION III
            
            p = sqrt(((w/cl)^2)-((w/cp(ii))^2));
            q = sqrt(((w/ct)^2)-((w/cp(ii))^2));
            val(ii) = ((tan(q*h))/q)+(4*(k(ii)^2)*p.*tan(p*h))/(((q.^2)-(k(ii)^2))^2);

        end
    end
end