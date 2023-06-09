function result = tensordot(A,B)
    result = zeros(3,3);
    for i=1:3
        for j=1:3
            for k = 1:3
                for l = 1:3
                    result(i,j)=result(i,j)+A(i,j,k,l)*B(k,l);
                end
            end
        end
    end
end