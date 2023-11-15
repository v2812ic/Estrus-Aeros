function [u, R] = solveSys(n_k, KG, Fext)

% The first node is embedded to the fuselage and its movement and rotation
% are restricted 

vR = [1, 2];
uR = [0, 0];

v_ndof = 1:(2*n_k + 1);
vL = setdiff(v_ndof, vR);

% Shorter matrices creation
KLL = KG(vL,vL);
KLR = KG(vL,vR);
KRL = KG(vR,vL);
KRR = KG(vR,vR);

% Computation of free and restricted forces
FextL = Fext(vL);
FextR = Fext(vR);

% System solver
uL = linsolve(KLL, FextL - KLR*uR');
R = KRR*uR' + KRL*uL - FextR;

% Computation of displacements vector
u = zeros(size(uR,1)+size(uL,1),1);
u(vL) = uL;
u(vR) = uR;


end