function result = ifourier3(expr,var,transvar)
    expr = ifourier(expr,var(1),transvar(1));
    expr = ifourier(expr,var(2),transvar(2));
    expr = ifourier(expr,var(3),transvar(3));
    result = expr;
end