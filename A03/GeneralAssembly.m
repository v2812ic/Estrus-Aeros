function [Fext, KG] = GeneralAssembly(n_k, Td, Kel, Fel)

Fext = zeros(2*(n_k+1), 1);
KG = zeros(2*(n_k + 1));

for i = 1 : n_k
    for j = 1 : 4 % 4 local dofs
        I = Td(i, j);
        Fext(I) = Fext(I) + Fel(j, i);
        for k = 1 : 4 
            J = Td(i, k);
            KG(I, J) = KG(I, J) + Kel(j, k, i);
        end
    end
end

end