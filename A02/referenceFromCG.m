function cg = referenceFromCG(n_dof,n_el,n_i,mat,Tmat,x,Tn,g)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_dof      Number of degrees of freedom
%                  n_el       Total number of elements
%                  n_i        Porblem's dimension
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
w_vec_1 = zeros(1,n_i);
w_vec_2 = zeros(1,n_i);
w_vec_T = zeros(1,n_i);
cg = zeros(size(x,1),size(x,2));
w_T = 0;

for i = 1 : n_el
    % Find the vector and length for each element
    le_vec = transpose(x(Tn(i,2),:)-x(Tn(i,1),:));
    le = norm(le_vec);

    % Calculate weight of the bar
    w = mat(Tmat(i),3)*mat(Tmat(i),2)*le*g;
    w_T = w_T + w;

    % Generate weight vector for each node and add them up
    w_vec_1 = w/2*x(Tn(i,1),:);
    w_vec_2 = w/2*x(Tn(i,2),:);
    w_vec_T = w_vec_T + w_vec_1 + w_vec_2;

end

w_vec_T = w_vec_T/w_T;

for i = 1 : size(cg,1)
    % Modify reference point
    cg(i,:) = x(i,:)-w_vec_T;
end

end