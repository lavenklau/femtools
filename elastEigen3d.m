C = sym('C',[3,3,3,3]);
xi = sym('xi',[3,1]);

if 0
% anisotroy
    Cmat = sym("C",[6,6]);
    assume(Cmat,'real');
    assume(xi,'real');
    for i = 1 : 6
        for j = i : 6
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

% symmetry
idtrans=[1, 6, 5; 
         6, 2, 4; 
         5, 4, 3];
for i = 1 : 3
    for j = 1 : 3
        for k = 1 : 3
            for l = 1 : 3
                C(i,j,l,k)=Cmat(idtrans(i,j),idtrans(k,l));
            end
        end
    end
end


A = sym(zeros(3));


for i = 1 : 3
    for k = 1 : 3
        for j = 1 : 3
            for l = 1 : 3
                A(i,k) = A(i,k) + C(i,j,k,l) * xi(j) * xi(l);
            end
        end
    end
end

% isotropic material
% sbv = [Cmat(1,1), lam+2*mu,     Cmat(1,2),lam,       Cmat(1,3),lam, ...
%            Cmat(2,1), lam,      Cmat(2,2),lam+2*mu,  Cmat(2,3),lam, ...
%            Cmat(3,1), lam,      Cmat(3,2),lam,       Cmat(3,3),lam+2*mu, ...
%            Cmat(4,4), mu,       Cmat(5,5),mu,        Cmat(6,6),mu];
% old = sbv(1:2:length(sbv));
% new = sbv(2:2:length(sbv));
% A = subs(A,old,new);
% for j = 4:6
%     for i = 1:6
%         A = subs(A,Cmat(i,j),0);
%     end
% end

Aj = adjoint(A);

[V,D]=eig(A);

Vsim = V * diag([xi(1),xi(1),xi(3)]);

% Vsim(:,2) = cross(Vsim(:,1),Vsim(:,3));

% Vsim(:,1) = cross(Vsim(:,2),Vsim(:,3));




