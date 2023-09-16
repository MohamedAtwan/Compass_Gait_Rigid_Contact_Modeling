function bm = modifyPolynomial(b, k)
bm = [];
for i = 1:k+1
    if i <= size(b,2)
        bm = [bm b(i)];
    else
        bm = [0 bm];
    end
end
