function [dc, omega_eq_2, beta_eq] = solve_c(t, c, G, beta, omega0, beta0, S0, b0, epsilon, alpha)
    omega_A = @(A) (omega0*sqrt(1 + 0.75*epsilon*A.^2));

    %beta eq
    int_beta_eq = @(A) ((A./(omega_A(A).^(1-alpha))).*exp((-G.*A.^2)/(2*c)));
    int_beta_eq_result = integral(int_beta_eq, 0, Inf);
    beta_eq = -beta0 + ((beta*G*sin(alpha*pi/2))/c)*int_beta_eq_result;

    % figure(1)
    % hold on
    % plot(t, beta_eq + beta0, '.r')
    % title('Beta equivalent')

    %omega eq
    first_int = @(A) (((omega_A(A)).^alpha).*A.*exp((-G.*A.^2)/(2*c)));
    first_int_result = integral(first_int, 0, Inf);
    second_int = @(A) (A.^3.*exp((-G.*A.^2)/(2*c)));
    second_int_result = integral(second_int, 0, Inf);
    omega_eq_2 = omega0^2 + ((beta*G*cos(alpha*pi/2))/c)*first_int_result + ((3*epsilon*(omega0^2)*G)/(4*c))*second_int_result;
    % 
    % figure(2)
    % hold on
    % plot(t, sqrt(omega_eq_2), '.k');
    % title('Omega equivalent')
    % 
    %eps
    Sw = evolutionary_power_spectrum(sqrt(omega_eq_2), t, S0, b0);
    
    %c(t)
    dc = -(beta0 + beta_eq)*c + pi*G*Sw/omega_eq_2;