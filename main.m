%% strucutre data
clc
clear

ndof = 3;

m = 1*ones(1,ndof); %mass vector
damping = 0.1*ones(1,ndof); %damping vector
k = 10*ones(1,ndof); %stiffness vector
ex = 0.5*ones(1,ndof);
ev = 0.2*ones(1,ndof);
T = 30;
b0 = 0.15;
S0 = 1;
alpha = 0.75;

[M, C, K] = get_mck(m, damping, k, ndof);

%% statistical_linearization
ntime = 200;
nfreq = 1000;
time = linspace(0, T, ntime);
omega_natural = sqrt(eig(inv(M)*K));
freq = linspace(0,max(omega_natural)*3,nfreq);

[var_displacement, var_velocity, conv, k_eq, c_eq] = statistical_linearization(m, damping, k, M, C, K, freq, time, ndof, S0, b0, ex, ev, alpha);

%% stochastic averaging
omega_eq_2 = var_velocity./var_displacement;
c = var_displacement;
dc(1,:) = gradient(c(1,:), time);
dc(2,:) = gradient(c(2,:), time);
dc(3,:) = gradient(c(3,:), time);
beta_eq = pi.*evolutionary_power_spectrum(sqrt(omega_eq_2), time, S0, b0)./(omega_eq_2.*c) - dc./c;

%% c(t)
beta0 = 0.07;
epsilon = 0.5;
tspan = time;

%analytical
ic = 0.00000001;
c_stored = zeros(ndof, numel(time));

for i=1:ndof
    beta_eq_dof = beta_eq(i,:);
    omega_eq_2_dof = omega_eq_2(i,:);

    beta_eq_dof(1) = beta_eq_dof(2);
    omega_eq_2_dof(1) = omega_eq_2_dof(2);

    G = sin(alpha*pi/2)./(omega_eq_2_dof.^(1-alpha));
    [t, c_stored(i,:)] = ode45(@(t, c_aux) solve_c_mdof(t, c_aux, G, beta_eq_dof, omega_eq_2_dof, S0, beta0, b0, time), tspan, ic);
end

plot(time, c_stored(1,:)', time, c(1,:)')
legend('analytical', 'sl')

%mcs
ns = 60;
var_x = 0; %displacement_variance_mcs(omega0, epsilon, beta0, alpha, S0, b0, tmax, ns);

%% plot
figure(2)
plot(t,c/G,t,var_x)
title("Oscillator non-stationary response variance E[x^2] = G^{-1}c(t)")

figure(3)
plot(t,sqrt(omega_eq_2))
title("\omega_{eq}(t)")

figure(4)
plot(t,beta_eq + beta0)
title("\beta_{eq}(t) + \beta_0")

%% survival probability
barrier = 0.35;

r_2 = zeros(numel(t)-1, 1);
time_domain = zeros(numel(t), 1);
time_domain(1) = t(1);
q = 0.07;

for i = 2:numel(time_domain)
    time_domain(i) = time_domain(i-1) + q*2*pi/sqrt((omega_eq_2(i-1)));
end

c_new = interp1(t, c, time_domain, 'pchip');
beta_eq_new = interp1(t, beta_eq, time_domain, 'pchip');
omega_eq_2_new = interp1(t, omega_eq_2, time_domain, 'pchip');
r_2(1) = 0;

for i = 2:numel(time_domain)
    tau = time_domain(i) - time_domain(i-1);
    r_2(i) = (c_new(i-1)/c_new(i))*(1 - (beta0 + beta_eq_new(i-1))*tau);
end

N = 5;
k = 1:1:N;
for i = 2:numel(r_2)
    A = (-G*barrier^2)/(2*c_new(i)*(1-r_2(i)));
    B = (-G*barrier^2)/(2*c_new(i-1)*(1-r_2(i)));
    D0 = (1-r_2(i))*exp(A)*(1 - exp(B));

    sum = 0;
    for n=1:N
        A = (r_2(i)^n)*(G^(2*n+2));
        B = (c_new(i-1)*c_new(i))^(n+1);
        C = ((1 - r_2(i))^(2*n+1))*prod((2*(1:n)).^2);
        upper_incomplete_gamma_i = igamma(n+1,(G*barrier^2)/(2*c_new(i)*(1-r_2(i))));
        upper_incomplete_gamma_i_minus = igamma(n+1,(G*barrier^2)/(2*c_new(i-1)*(1-r_2(i))));
        complete_gamma = gamma(n+1);
        Ln = ((4^n)*((1 - r_2(i))^(2*n+2))*(c_new(i-1)^(n+1))*(c_new(i)^(n+1))*(G^(-2*n-2))*upper_incomplete_gamma_i)*(complete_gamma - upper_incomplete_gamma_i_minus);
        Dn = (A/(B*C))*Ln;
        sum = sum + Dn;
    end
    Q(i,1) = D0 + sum;
    H(i,1) = 1 - exp((-G*barrier^2)/(2*c_new(i-1)));
    F(i,1) = Q(i,1)/H(i,1);
end

for i = 1:numel(time_domain)
    aux = 1;
    for M = 1:i
        aux = aux*(1 - F(M));
    end
    Pb(i, 1) = aux;
end

figure(5)
plot(time_domain, Pb)
xlim([0 15])