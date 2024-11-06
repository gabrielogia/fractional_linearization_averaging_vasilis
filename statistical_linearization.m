function [var_displacement, var_velocity, conv, k_eq_time, c_eq_time] = statistical_linearization(m, c, k, M, C, K, freq, time, ndof, S0, b0, ex, ev, alpha)
    Mt = M;
    tol = 1e-6;
    maxiter = 30;
    ntime = numel(time);
    nfreq = numel(freq);
    q = alpha;

    var_displacement = zeros(ndof,ntime);
    var_velocity = zeros(ndof,ntime);
    k_eq_time = zeros(ndof,ntime); 
    c_eq_time = zeros(ndof,ntime);
    keq = zeros(ndof, 1);
    ceq = zeros(ndof, 1);

    for i=1:ntime %loop in time
        t = time(i);
        [Ceq, Keq] = get_equivalent_ck(ceq, keq, ndof);

        sx2 = zeros(ndof, 1);
        sv2 = zeros(ndof, 1);
    
        dkeq = 1000*ones(ndof,1);
        dceq = 1000*ones(ndof,1);
        dkeq_max = max(dkeq);
        dceq_max = max(dceq);

        itera = 0;
        while (dkeq_max > tol || dceq_max > tol) && (itera < maxiter)
            Ct = C + Ceq;
            Kt = K + Keq;
    
            H_ps = zeros(ndof, nfreq);
            H_ps_freq = zeros(ndof, nfreq);
          
            for j=1:nfreq %loop in frequency
                f = freq(j);
                H = get_H(f,Mt,Ct,Kt, q);
                
                aux = zeros(ndof,1);
                for p = 1:ndof
                    ps = evolutionary_power_spectrum(f, t, S0, b0);
                    aux = aux + (2*ps*(m(p).^2)*abs(H(1:ndof, p)).^2);
                end

                H_ps(1:ndof, j) = aux;
                H_ps_freq(1:ndof, j) = (f.^2)*aux;
            end
    
            ceq_1 = ceq;
            keq_1 = keq;
    
            sx2 = zeros(ndof, 1);
            sv2 = zeros(ndof, 1);
            ceq = zeros(ndof, 1);
            keq = zeros(ndof, 1);
         
            for l=1:ndof
                Ex = trapz(freq, H_ps(l,:));
                Exd = trapz(freq, H_ps_freq(l,:));
             
                sx2(l) = Ex;
                sv2(l) = Exd;
                ceq(l) = 3*ev(l)*c(l)*Exd;
                keq(l) = 3*ex(l)*k(l)*Ex;
            end
            
            [Ceq, Keq] = get_equivalent_ck(ceq, keq, ndof);
     
            dceq = abs(ceq - ceq_1)./ceq_1;
            dkeq = abs(keq - keq_1)./keq_1;

            conv(i).ceq(:,itera+1) = abs(ceq - ceq_1);
            conv(i).keq(:,itera+1) = abs(keq - keq_1);

            dceq_max = max(dceq);
            dkeq_max = max(dkeq);
            
            itera = itera + 1;
        end

        var_displacement(:,i) = sx2;
        var_velocity(:,i) = sv2;
        k_eq_time(:,i) = keq;
        c_eq_time(:,i) = ceq;
    end