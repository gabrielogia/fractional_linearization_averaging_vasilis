function [dc, omega_eq_2, beta_eq] = solve_c(t, c, G, beta, omega0, beta0, S0, b0, epsilon, alpha)
    omega_A = @(A) (omega0*sqrt(1 + 0.75*epsilon*A.^2));

    %beta eq
    int_beta_eq = @(A) ((A./(omega_A(A).^(1-alpha))).*exp((-G.*A.^2)/(2*c)));
    int_beta_eq_result = integral(int_beta_eq, 0, Inf);
    beta_eq = -beta0 + ((beta*G*sin(alpha*pi/2))/c)*int_beta_eq_result;

    %omega eq
    first_int = @(A) (((omega_A(A)).^alpha).*A.*exp((-G.*A.^2)/(2*c)));
    first_int_result = integral(first_int, 0, Inf);
    second_int = @(A) (A.^3.*exp((-G.*A.^2)/(2*c)));
    second_int_result = integral(second_int, 0, Inf);
    omega_eq_2 = omega0^2 + ((beta*G*cos(alpha*pi/2))/c)*first_int_result + ((3*epsilon*(omega0^2)*G)/(4*c))*second_int_result;
    
    %eps
    aux = omega_eq_2/((5*pi)^2);
    Sw = S0 * aux * exp(-b0*t)*t^2 * exp(-aux*t);
    
    %c(t)
    dc = -(beta0 + beta_eq)*c + pi*G*Sw/omega_eq_2;