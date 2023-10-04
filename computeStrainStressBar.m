function [eps,sig] = computeStrainStressBar(n_d,n_el,u,Td,x,Tn,mat,Tmat,Delta_T)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_d        Problem's dimensions
%                  n_el       Total number of elements
%   - u     Global displacement vector [n_dof x 1]
%            u(I) - Total displacement on global DOF I
%   - Td    DOFs connectivities table [n_el x n_el_dof]
%            Td(e,i) - DOF i associated to element e
%   - x     Nodal coordinates matrix [n x n_d]
%            x(a,i) - Coordinates of node a in the i dimension
%   - Tn    Nodal connectivities table [n_el x n_nod]
%            Tn(e,a) - Nodal number associated to node a of element e
%   - mat   Material properties table [Nmat x NpropertiesXmat]
%            mat(m,1) - Young modulus of material m
%            mat(m,2) - Section area of material m
%   - Tmat  Material connectivities table [n_el]
%            Tmat(e) - Material index of element e
%   - Delta_T Temperature increment causing the deformation
%--------------------------------------------------------------------------
% It must provide as output:
%   - eps   Strain vector [n_el x 1]
%            eps(e) - Strain of bar e
%   - sig   Stress vector [n_el x 1]
%            sig(e) - Stress of bar e
%--------------------------------------------------------------------------

% Initialize
ue = zeros(n_d*2,1);
Re = zeros(2,n_d*2);
eps = zeros(n_el,1);
sig = zeros(n_el,1);

% Iterate over all elements
for i = 1 : n_el
    % Find length vector components
    le_vec = x(Tn(i,2),:)-x(Tn(i,1),:);
    le = norm(le_vec);
    
    % Rotation matrix
    Re(1,1:n_d) = le_vec;
    Re(2,n_d+1:size(Re,2)) = le_vec;

    % Apply rotation
    ue = u(Td(i,:));
    uep = (Re*ue)/le;

    % Find epsilon
    aux = [-1 1];
    eps_e = (aux*uep)/le;

    % Find sigma
    Ee = mat(Tmat(i),1);
    sig_e = Ee*eps_e;

    % Store values
    eps(i) = eps_e;
    sig(i) = sig_e;

end

end