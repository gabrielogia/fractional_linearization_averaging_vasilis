function Sw = evolutionary_power_spectrum(freq, t, S0, b0)
    aux = (freq./(5*pi)).^2;
    Sw = S0.*aux.*(exp(-b0.*t).*t.^2).*exp(-aux.*t);