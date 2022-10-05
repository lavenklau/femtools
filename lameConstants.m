% compute lame constants from young's modulus and poisson ratio
%
%  Author : lavenklau
%  mail   : zhd9702@mail.ustc.edu.cn
function [lam, mu] = lameConstants(E, v)
    lam = v * E / ((1 + v)*(1 - 2 * v));
    mu  = E / (2*(1 + v));
end
