function [pu,pt,Fy,Mz] = findProperties(n_i,n_k,Td,Tnod,Kel,u)

ue = zeros(2*n_i);
Fy = zeros(n_k,2);
Mz = zeros(n_k,2);
Coeffs_matrix = zeros(4,4);
pu = zeros(n_k,4);
pt = zeros(n_k,3);
abcd = zeros(4);

for i = 1 : n_k
    le = Tnod(i,2) - Tnod(i,1);
    ue = u(Td(i,:));

    % Auxiliar matrix assembly
    Coeffs_matrix = [[2,le,-2,le];
                     [-3*le,-2*le^2,3*le,-le^2];
                     [0,le^3,0,0];
                     [le^3,0,0,0]];
    Coeffs_matrix = Coeffs_matrix/(le^3);

    % Internal loads
    Loads_e = Kel(:,:,i)*ue;

    % Shear force
    Fy(i,1) = -Loads_e(1);
    Fy(i,2) = Loads_e(3);
    
    % Bending moment
    Mz(i,1) = -Loads_e(2);
    Mz(i,2) = Loads_e(4);

    % Third order polynomial coefficients
    abcd = Coeffs_matrix*ue;

    % Displacement & rotation coefficients, respectively
    pu(i,:) = abcd;
    pt(i,:) = [3*abcd(1), 2*abcd(2), abcd(3)];
end

end