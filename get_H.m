function H = get_H(freq, Mt, Ct, Kt, q)
    H = inv(-(freq.^2) * Mt + (1j * freq)^q * Ct + Kt);
end