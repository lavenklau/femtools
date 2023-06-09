function result = Tetbern3d(p, i, j, k, xi, eta, zet)
    c = nchoosek(p,i) * nchoosek(p-i,j) * nchoosek(p-i-j,k);
    x = xi.^i .* eta.^j .* zet.^k .* (1-xi-eta-zet).^(p-i-j-k);
    result = c*x;
end