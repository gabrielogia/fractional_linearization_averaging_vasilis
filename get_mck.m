function [M, C, K] = get_mck(m, c, k, ndof)
    M = zeros(ndof, ndof);
    C = zeros(ndof, ndof);
    K = zeros(ndof, ndof);
    
    for i=1:ndof
        C(i,i) = c(i);
        K(i,i) = k(i);

        if i < ndof
            C(i,i+1) = -c(i+1);
            K(i,i+1) = -k(i+1);
        end

        for j=1:i
            M(i,j) = m(i);
        end
    end