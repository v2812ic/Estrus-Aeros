function pu = findDisplacementCoefficients(n_k,Tnod)

ue = zeros()

for i = 1 : n_k
    le = Tnod(i, 2) - Tnod(i, 1);
    ue(:,i)
end