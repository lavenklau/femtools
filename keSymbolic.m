% Param
%  -input :
%    @uselame : whether to represent the symbolic matrix via lame 
%               constants (lam, mu), if not, use youngs modulus and poisson
%               ratio (E, v).
%    @unitElement : whether the element has unit length along all axis.
%
%  -output:
%    @kesymb : symbolic element stiffness matrix representation
%    @kelam,kemu: matrix for each lame constant, the stiffness matrix can
%              be represent by KE = lam * kelam + mu * kemu
function [kesymb, kelam, kemu] = keSymbolic(uselame,unitElement)

if ~unitElement
    len = sym("l%d%d",[3,1]); 
else
    len = sym([1;1;1]);
end

x = sym('x%d%d',[3,1]);

assume(len,'positive');
assume(x,'real');

i0 = 1-x(1)/len(1);
i1 = x(1)/len(1);
j0 = 1-x(2)/len(2);
j1 = x(2)/len(2);
k0 = 1-x(3)/len(3);
k1 = x(3)/len(3);

N = [i0 * j0 * k0; ...
     i1 * j0 * k0; ...
     i0 * j1 * k0; ...
     i1 * j1 * k0; ...
     i0 * j0 * k1; ...
     i1 * j0 * k1; ...
     i0 * j1 * k1; ...
     i1 * j1 * k1;];
dN = [diff(N,'x11'), diff(N,'x21'), diff(N,'x31')];

B = [ ...
   bi(1,dN), bi(2,dN), bi(3,dN), bi(4,dN), ...
   bi(5,dN), bi(6,dN), bi(7,dN), bi(8,dN)
];

% expression in young's modulus and poisson ratio
E = sym('E');
v = sym('v');
% expression in lame constants
lam = sym('lam');
mu  = sym('mu'); 

if ~uselame
    k = E * (1 - v) / ((1 + v) * (1 - 2 * v));
    mud = v / (1 - v);
    mut = (1 - 2 * v) / (2 * (1 - v));
    D = [   [1, mud, mud, 0, 0, 0]; ...
            [mud, 1, mud, 0, 0, 0]; ...
            [mud, mud, 1, 0, 0, 0]; ...
            [0, 0, 0, mut, 0, 0]; ...
            [0, 0, 0, 0, mut, 0]; ...
            [0, 0, 0, 0, 0, mut]] * k;
else
    D = lam * ...
      [ [1, 1, 1, 0, 0, 0]; ...
        [1, 1, 1, 0, 0, 0]; ...
        [1, 1, 1, 0, 0, 0]; ...
        [0, 0, 0, 0, 0, 0]; ...
        [0, 0, 0, 0, 0, 0]; ...
        [0, 0, 0, 0, 0, 0]  ...
      ] + ...
        mu * [  ...
        [2, 0, 0, 0, 0, 0]; ...
        [0, 2, 0, 0, 0, 0]; ...
        [0, 0, 2, 0, 0, 0]; ...
        [0, 0, 0, 1, 0, 0]; ...
        [0, 0, 0, 0, 1, 0]; ...
        [0, 0, 0, 0, 0, 1]  ...
      ];
end

ke = int(B'*D*B,'x11',0,len(1));
ke = int(ke,'x21',0,len(2));
ke = int(ke, 'x31',0,len(3));
kesymb = ke;
if ~uselame
    % kev = eval(subs(ke,{'E','v','l11','l21','l31'},[1e6,0.3,1,1,1]));
    kelam = [];
    kemu  = [];
else
    kelam = simplify(subs(ke,'mu',0)/lam);
    kemu  = simplify(subs(ke,'lam',0)/mu);
    %latx = latex(ke);
end

end

function result = bi(i,dN)
    result = [  ...
        dN(i,1), 0, 0;  ...
        0, dN(i,2), 0;  ...
        0, 0, dN(i,3);  ...
        0, dN(i,3), dN(i,2);  ...
        dN(i,3), 0, dN(i,1);  ...
        dN(i,2), dN(i,1), 0;  ...
    ];
end