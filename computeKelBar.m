function Kel = computeKelBar(n_k, Tnod, mat, Tmat)

Kel = zeros(4, 4, n_k); % Dimensions are set to 4x4

for i = 1 : n_k
    le = Tnod(i, 2) - Tnod(i, 1);
    f = mat(Tmat(i), 3)*mat(Tmat(i), 1)/(le^3);
    Kel(:, :, i) = f*[[12, 6*le, -12, 6*le];[6*le, 4*le^2, -6*le, 2*le^2];[-12, -6*le, 12, -6*le];[6*le, 2*le^2, -6*le, 4*le^2]];
end
end