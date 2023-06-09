%% Solve 1-D poisson problem with non-homogeneous boundary condition
% 
% $$\Delta u = - f, \quad x\in(0,1)$$
% $f = x(1-x)sin(x), u(0)=1, u(1)=2$ 
% 
% 

h = 0.001;

x = 0:h:1;

nv = length(x);
ne = length(x) - 1;

x0 = 1; x1 = 2;
ir = [];
jc = [];
jv = [];
b = zeros(n,1);
for e = 1 : ne
    u1 = e;
    u2 = e + 1;
    g1 = - 1 / h;
    g2 = 1 / h;
    % one point quadrature position
    p = (e -0.5) * h;
    b(u1) = 0.5 / h * p * (1-p) * sin(p); % 
    b(u2) = 0.5 / h * p * (1-p) * sin(p);
    if e == 1 
        g1 = 0;
    elseif e == ne
        g2 = 0;
    end
    ir = [ir ; e ; e ; e+1; e+1];
    jc = [jc ; e ; e+1 ; e ; e+1];
    jv = [jv ; g1 * g1 ; g1 * g2 ; g2 * g1 ; g2 * g2];
end



b(1) = 0;
b(n) = 0;
b(2) = - (1 / h * x0 * (-1 / h));
b(n-1) = - (- 1 / h * x1 * (1 / h));

K = sparse(ir, jc, jv);
K(1,1)=1;
K(nv,nv)=1;

u = K\b;
u(1) = u(1) + x0;
u(nv)= u(nv) + x1;

