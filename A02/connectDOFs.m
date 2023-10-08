function Td = connectDOFs(n_el,n_nod,n_i,n_el_dof,Tn)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_el     Total number of elements
%                  n_nod    Number of nodes per element
%                  n_i      Number of DOFs per node
%   - Tn    Nodal connectivities table [n_el x n_nod]
%            Tn(e,a) - Nodal number associated to node a of element e
%--------------------------------------------------------------------------
% It must provide as output:
%   - Td    DOFs connectivities table [n_el x n_el_dof]
%            Td(e,i) - DOF i associated to element e
%--------------------------------------------------------------------------

Td = zeros(n_el,n_el_dof);
    
for i = 1 : n_nod
    for j = 1 : n_i
        Td(:,(i - 1)*n_i + j) = Tn(:,i)*n_i-(n_i-j);
    end
end

end