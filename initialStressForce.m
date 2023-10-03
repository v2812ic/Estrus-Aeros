function Fint_0 = initialStressForce(n_dof,n_el,n_i,x,Tn,mat,Tmat,Delta_T);
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_i         Number of DOFs per node
%                  n_dof       Total number of DOFs
%                  n_el       Total number of elements
%   - x     Nodal coordinates matrix [n x n_d]
%            x(a,i) - Coordinates of node a in the i dimension
%   - Tn    Nodal connectivities table [n_el x n_nod]
%            Tn(e,a) - Nodal number associated to node a of element e
%   - mat   Material properties table [Nmat x NpropertiesXmat]
%            mat(m,1) - Young modulus of material m
%            mat(m,2) - Section area of material m
%   - Tmat  Material connectivities table [n_el]
%            Tmat(e) - Material index of element e
%   - Delta_T        Temperature increment causing the deformation
%--------------------------------------------------------------------------
% It must provide as output:
%   - Fint_0  Initial stress force vector [n_dof x 1]
%             Fint_0(I) - Internal force acting on DOF I
%--------------------------------------------------------------------------

Fint_0 = zeros(n_dof,1);

% Create the internal forces due to the thermal elongation
for i = 1 : n_el
    F0_vec = zeros(n_dof,1);
    % Find the vectors and inverse vector for each element
    le_vec = transpose(x(Tn(i,2),:)-x(Tn(i,1),:));
    le = norm(le_vec);

    % Force from every element
    F0 = mat(Tmat(i),2)*mat(Tmat(i),1)*mat(Tmat(i),3)*Delta_T;
    F0_vec(Tn(i,2)*n_i-n_i+1:Tn(i,2)*n_i,1) = F0*le_vec/le;
    F0_vec(Tn(i,1)*n_i-n_i+1:Tn(i,1)*n_i,1) = -F0*le_vec/le;

    % Accumulated force
    Fint_0 = Fint_0 + F0_vec;
end

end