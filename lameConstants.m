function [lam, mu] = lameConstants(E, v)
    lam = v * E / ((1 + v)*(1 - 2 * v));
    mu  = E / (2*(1 + v));
end