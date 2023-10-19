function [L,D,T] = computeLD(n,Fw,Fext,x,Aeronodes,Tnodes)
%--------------------------------------------------------------------------
% The function takes as inputs:
%  - Dimensions: n  Number of nodes
%  - x              Nodal coordinates matrix
%  - Fw, Fext       Mass and external forces
%  - Aeronodes      Nodes where aerodynamic forces are applied
%  - Tnodes         Nodes where thrust is applied
%--------------------------------------------------------------------------
% It must provide as output:
%  - L              Total Lift
%  - D              Total Drag
%  - T              Total Thrust
% n_d wasn't included since a particular case is being considered 
% in this code (3D structure).
% Also this only works if the structure and the forces are symmetric
% Fint_0 is not considered since it doesn't generate net momentum
%--------------------------------------------------------------------------

A = zeros(3); 
B = zeros(3, 1);

%Horizontal and vertical force equilibrium
A(1, 1) = 1; %V
A(2, 2) = 1; %H
A(2, 3) = 1; %H
for i = 1:n
    B(1) = B(1) - Fext(3*i) - Fw(3*i);%V
    B(2) = B(2) - Fext(3*i - 2) - Fw(3*i - 2);%H
end

%y-axis momentum equilibrium
A(3, 1) = sum(x(Aeronodes, 1))/size(Aeronodes, 2);
A(3, 2) = -sum(x(Aeronodes, 3))/size(Aeronodes, 2);
A(3, 3) = -sum(x(Tnodes, 3))/size(Tnodes, 2); 
for i = 1:n
    B(3) = B(3) + x(i, 3)*(Fext(3*i-2)+Fw(3*i-2)) - x(i, 1)*(Fext(3*i)+Fw(3*i));
end

% Solve system and assign values
C = A\B;
L = C(1);
D = C(2);
T = C(3);
end