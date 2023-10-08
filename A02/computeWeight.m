function Fw = computeWeight(n_dof,n_el,mat,Tmat,x,Tn,g)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_dof      Number of degrees of freedom
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
%   - g     Gravity
%--------------------------------------------------------------------------
% It must provide as output:
%   - Fw   Weight force from each bar to each node
%--------------------------------------------------------------------------

Fw = zeros(n_dof,1);

for i = 1 : n_el
    % Find the vector and length for each element
    le_vec = transpose(x(Tn(i,2),:)-x(Tn(i,1),:));
    le = norm(le_vec);

    % Calculate weight
    w = mat(Tmat(i),3)*mat(Tmat(i),2)*le*g;

    % Imposing the force over each node of the bar
    Fw(Tn(i,1)*3) = Fw(Tn(i,1)*3) - w/2;
    Fw(Tn(i,2)*3) = Fw(Tn(i,2)*3) - w/2;

end

end