function [Fext, KG] = GeneralAssembly(n_k, Td, Kel, Fel, Me, g, nelk, L1, L)

Fext = zeros(2*(n_k + 1), 1);
KG = zeros(2*(n_k + 1));
dofmot = L1/L*nelk*2+1;

for i = 1 : n_k
    for j = 1 : 4 % 4 local dofs
        % meter la fuerza en el grado de libertad
        I = Td(i, j);
        % fuerza del motor, considerando que solo puede caer sobre un nodo,
        % no sobre una viga general
        if I == dofmot
            Fext(I) = Fext(I) - Me*g/2;   
            fprintf('hola concha');
        end
        Fext(I) = Fext(I) + Fel(j, i);
        for k = 1 : 4 
            J = Td(i, k);
            KG(I, J) = KG(I, J) + Kel(j, k, i);
        end
    end
end

end