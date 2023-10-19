function cg = referenceFromCG(n,n_el,n_i,m_p,mat,Tmat,x,Tn)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_dof      Number of degrees of freedom
%                  n_el       Total number of elements
%                  n_i        Problem's dimension
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
%   - cg    Nodal coordinates from the center of gravity
%--------------------------------------------------------------------------

m_vec_1 = zeros(1,n_i);
m_vec_2 = zeros(1,n_i);
m_vec_T = zeros(1,n_i);
cg = zeros(size(x,1),size(x,2));
m_T = 0;

for i = 1 : n_el
    % Find the vector and length for each element
    le_vec = transpose(x(Tn(i,2),:)-x(Tn(i,1),:));
    le = norm(le_vec);

    % Calculate the meight of the bar
    m = mat(Tmat(i),3)*mat(Tmat(i),2)*le;
    m_T = m_T + m;

    % Generate the meight vector for each node and add them up
    m_vec_1 = m/2*x(Tn(i,1),:);
    m_vec_2 = m/2*x(Tn(i,2),:);
    m_vec_T = m_vec_T + m_vec_1 + m_vec_2;

end

for i = 1 : n
    % Add passenger's mass
    m_vec_T = m_vec_T + m_p(i)*x(i,:);

end

% Generate final vector
m_T = m_T + sum(m_p);
m_vec_T = m_vec_T/m_T;

for i = 1 : size(cg,1)
    % Modify reference point
    cg(i,:) = x(i,:)-m_vec_T;
end

end