function [u,R] = solveSys(vL,vR,uR,KG,Fext)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - vL      Free degree of freedom vector
%   - vR      Prescribed degree of freedom vector
%   - uR      Prescribed displacement vector
%   - KG      Global stiffness matrix [n_dof x n_dof]
%              KG(I,J) - Term in (I,J) position of global stiffness matrix
%   - Fext    Global force vector [n_dof x 1]
%              Fext(I) - Total external force acting on DOF I
%--------------------------------------------------------------------------
% It must provide as output:
%   - u       Global displacement vector [n_dof x 1]
%              u(I) - Total displacement on global DOF I
%   - R       Global reactions vector [n_dof x 1]
%              R(I) - Total reaction acting on global DOF I
%--------------------------------------------------------------------------

% Shorter matrices creation
KLL = KG(vL,vL);
KLR = KG(vL,vR);
KRL = KG(vR,vL);
KRR = KG(vR,vR);

% Computation of free and restricted forces
FextL = Fext(vL);
FextR = Fext(vR);

% System solver
uL = inv(KLL)*(FextL - KLR*uR);
R = KRR*uR + KRL*uL - FextR;

% Computation of displacements vector
u = zeros(size(uR,1)+size(uL,1),1);
u(vL) = uL;
u(vR) = uR;

end