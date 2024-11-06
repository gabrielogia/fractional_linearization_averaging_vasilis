function dc = solve_c_mdof(t, c, G, beta, omega0, S0, beta0, b0, time)
    %beta eq
    beta_eq = interp1(time, beta, t);

    %omega eq
    omega_eq_2 = interp1(time, omega0, t);

    %G in time
    g = interp1(time, G, t);
   
    %eps
    Sw = evolutionary_power_spectrum(sqrt(omega_eq_2), t, S0, b0);
    
    %c(t)
    dc = -(beta0 + beta_eq)*c + pi*g*Sw/omega_eq_2;