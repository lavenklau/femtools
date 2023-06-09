%% L maps coefficients of bernstein to nodal values, its inverse maps nodal 
%  values to bernstein coefficients
L=[];
r = 0;
c = 0;
for zet = 0 : 0.5 : 1
    for eta = 0 : 0.5 : 1 - zet
        for xi = 0 : 0.5 : 1-eta-zet
            r = r + 1;
            c = 0;
            for i = 0 : 1 : 2
                for j = 0 : 1 : 2 - i
                    for k = 0 : 1 : 2-j-i
                        c = c + 1;
                        L(r,c) = Tetbern3d(2,i,j,k,xi,eta,zet);
                    end
                end
            end
        end
    end
end

%% verity tangent moduli
F = [1, 0.1, 0; 0.1, 1, 0.1; 0, 0.1, 1.1];

PK2_old = evalPK2(F);
PK1_old = F * evalPK2(F);
E_old = (F'*F-eye(3))/2;
phi_old = evalPhi(F);
PK1_num = [];
ctan = evalCtan(F);

for i = 1:3
    for j = 1:3
        dF = zeros(3,3);
        dF(i,j)=1e-7;
        F1_new = F + dF;
        E_new = (F1_new'*F1_new-eye(3))/2;
        dE = (dF'*F+F'*dF)/2;
        PK1_new = (F + dF) * evalPK2(F + dF);
        PK2_new = evalPK2(F+dF);
        dPK2 = PK2_new - PK2_old;
        dPK1 = PK1_new - PK1_old;
        K_num(1:3,1:3,i,j) = (PK1_new - PK1_old) / 1e-7;
        PK1_num(i,j) = (evalPhi(F+dF)-phi_old)/1e-7;
    end
end

KdF   = zeros(3,3);
KdF58 = zeros(3,3);
KdF56 = zeros(3,3);
KdF57 = zeros(3,3);
KdF55 = zeros(3,3);
KdF55_2 = zeros(3,3);
ctandE = tensordot(ctan,dE);
for i=1:3
    for j=1:3
        for k=1:3
            KdF58(i,j) = KdF58(i,j) + dF(i, k) * PK2_old(k, j);
            KdF56(i,j) = KdF56(i,j) + dF(i,k)*PK2_old(k,j);
            KdF57(i,j) = KdF57(i,j) + dF(i,k)*PK2_old(k,j);
            for l=1:3
                for p=1:3
                    for q=1:3
                        KdF58(i,j)=KdF58(i,j) + ...
                        ctan(p,j,l,q)*F(k,q)*F(i,p)*dF(k,l);

                        KdF56(i,j)=KdF56(i,j) + ...
                        ctan(p,j,k,l)*(dF(q,k)*F(q,l)+F(q,k)*dF(q,l))* F(i,p)/2 ;

                        KdF57(i,j)=KdF57(i,j) + ...
                        0.5*(ctan(p,j,l,q)*F(k,q)*F(i,p)*dF(k,l) + ...
                        ctan(p,j,q,l)*F(k,q)*F(i,p)*dF(k,l));
                    end
                end
            end
        end
    end
end

K = evalPseudoModuli(F);


%%
fprintf('F*ctandE+dF*PK2_old = ');
disp(F*ctandE+dF*PK2_old);
tmp = zeros(3,3);
for i =1:3
    for j=1:3
        for k=1:3
            tmp(i,j) = tmp(i,j) + F(i,k)*ctandE(k,j) + dF(i,k)*PK2_old(k,j);
        end
    end
end
fprintf("tmp = ");
disp(tmp);

function result = evalPhi(F)
    C = F' * F;
    I_c = trace(C);
    J = det(F);logJ=log(J);
    lambda = 1e3;
    mu = 1e3;
    result = mu / 2 * (I_c - 3) - mu * log(J) + lambda / 2 * logJ^2;
end

function result = evalPK2(F)
    C = F' * F;
    Cinv = inv(C);
    mu = 1e3;
    lambda = 1e3;
    J = det(F);
    logJ = log(J);
    S = mu * (eye(3) - Cinv) + lambda * logJ * Cinv;
    result = S;
end

function result = evalCtan(F)
    C = F' * F;
    Cinv = inv(C);
    mu = 1e3;
    lambda = 1e3;
    J = det(F);
    logJ = log(J);
    S = mu * (eye(3) - Cinv) + lambda * logJ * Cinv;
    Ctan = [];
    I4 = [];
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    I4(i,j,k,l) = Cinv(i,k) * Cinv(j,l);
                end
            end
        end
    end

    I4 = (permute(I4,[1,2,3,4]) + permute(I4,[1,2,4,3]) + ...
        permute(I4,[2,1,3,4]) + permute(I4,[2,1,4,3]))/4;
    I4 = (permute(I4,[3,4,1,2]) + permute(I4,[3,4,1,2])) / 2;
    
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    Ctan(i,j,k,l) = lambda * Cinv(i,j) * Cinv(k,l) + ...
                    2 * (mu - lambda * logJ) * I4(i,j,k,l);
                end
            end
        end
    end

    result = Ctan;
end

function result = evalPseudoModuli(F)
    ctan = evalCtan(F);
    pk2 = evalPK2(F);
    result = zeros(3,3,3,3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    for p = 1:3
                        for q = 1:3
                            result(i,j,k,l) = result(i,j,k,l) + ...
                            ctan(p,j,l,q) * F(k,q) * ...
                            F(i,p) + pk2(l,j) * (i==k);
                        end
                    end
                end
            end
        end
    end
end


