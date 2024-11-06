function [Ceq, Keq]=get_equivalent_ck(ceq, keq, ndof)
    Ceq = zeros(ndof, ndof);
    Keq = zeros(ndof, ndof);

    for i=1:ndof
        Ceq(i, i) = ceq(i);
        Keq(i, i) = keq(i);
    end
