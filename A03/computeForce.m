function Fel = computeForce(n_k, Tnod, eqmass1, eqmass2, eqlift1, eqlift2, lp, L1, g)

Fel = zeros(4, n_k); % Dimensions are set to 4 x n_k
syms x l
eqlift1 = subs(eqlift1, l, lp);
eqlift2 = subs(eqlift2, l, lp);

for i = 1 : n_k
    x2 = Tnod(i, 2);
    x1 = Tnod(i, 1);
    le = x2 - x1;
    
    liftbar = (int(eqlift1, x, [min(x1, L1), min(x2, L1)]) + int(eqlift2, x, [max(x1, L1), max(x2, L1)]))/le;
    massbar = (int(eqmass1, x, [min(x1, L1), min(x2, L1)]) + int(eqmass2, x, [max(x1, L1), max(x2, L1)]))/le;
    
    qbar = double(liftbar - massbar*g);
    
    Fel(:, i) = qbar*le/2 * [1, le/6, 1, -le/6]; 
    
end