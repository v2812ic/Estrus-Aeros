function Kel = computeKelBar(n_d,n_el,n_el_dof,x,Tn,mat,Tmat)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_d        Problem's dimensions
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
%--------------------------------------------------------------------------
% It must provide as output:
%   - Kel   Elemental stiffness matrices [n_el_dof x n_el_dof x n_el]
%            Kel(i,j,e) - Term in (i,j) position of stiffness matrix for element e
%--------------------------------------------------------------------------

Kel = zeros(n_el_dof,n_el_dof,n_el); 
Re = zeros(2,n_d*2); 
Ke_p = [[1, -1]; [-1, 1]];

for i = 1 : n_el
    Ae = mat(Tmat(i), 2);
    Ee = mat(Tmat(i), 1);
    le_vec = x(Tn(i,2),:)-x(Tn(i,1),:);
    le = norm(le_vec);
    
    Re(1,1:n_d) = le_vec;
    Re(2,n_d+1:size(Re,2)) = le_vec;
    
    Kel(:,:,i) = Ae*Ee/(le^3)*transpose(Re)*Ke_p*Re;
end

end