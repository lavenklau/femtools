function KH = isosurfKH(ph, x)
    assume(x,'real');
    dx = gradient(ph,x);
    dim = length(x);
    ddx = sym(zeros(dim));
    for i = 1 : length(dx)
        ddx(1:dim,i) = gradient(dx(i),x);
    end
    n = dx / sqrt(dx'*dx);
    dx2 = dx'*dx;
    KH = simplify(divergence(n,x));
end





