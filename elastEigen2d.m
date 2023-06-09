C = sym('C',[2,2,2,2]);
xi = sym('xi',[2,1]);

if 0
% anisotroy
    Cmat = sym("C",[3,3]);
    assume(Cmat,'real');
    assume(xi,'real');
    for i = 1 : 3
        for j = i : 3
            Cmat(j,i) = Cmat(i,j);
        end
    end
else
% isotropy
    lam = sym('lambda');
    mu  = sym('mu');
    assume(lam,'real');
    assume(mu,'real');
    Cmat = lam * ...
          [ [1, 1, 0,]; ...
            [1, 1, 0,]; ...
            [0, 0, 0,]; ...
          ] + ...
            mu * [  ...
            [2, 0, 0,]; ...
            [0, 2, 0,]; ...
            [0, 0, 1,]; ...
          ];
end

% symmetry
idtrans=[1, 3; 
         3, 2];

for i = 1 : 2
    for j = 1 : 2
        for k = 1 : 2
            for l = 1 : 2
                C(i,j,l,k)=Cmat(idtrans(i,j),idtrans(k,l));
            end
        end
    end
end


A = sym(zeros(2));


for i = 1 : 2
    for k = 1 : 2
        for j = 1 : 2
            for l = 1 : 2
                A(i,k) = A(i,k) + C(i,j,k,l) * xi(j) * xi(l);
            end
        end
    end
end

Aj = adjoint(A);

[V,D]=eig(A);

Vsim = V * diag([xi(1),xi(2)]);

% Vsim(:,2) = cross(Vsim(:,1),Vsim(:,3));

% Vsim(:,1) = cross(Vsim(:,2),Vsim(:,3));