function [sig_cr] = analyseBuckling(n_el,x,Tn,mat,Tmat,sig)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_el       Total number of elements
%   - x     Nodal coordinates matrix [n x n_d]
%            x(a,i) - Coordinates of node a in the i dimension
%   - Tn    Nodal connectivities table [n_el x n_nod]
%            Tn(e,a) - Nodal number associated to node a of element e
%   - mat   Material properties table [Nmat x NpropertiesXmat]
%            mat(m,1) - Young modulus of material m
%            mat(m,2) - Section area of material m
%   - Tmat  Material connectivities table [n_el]
%            Tmat(e) - Material index of element e
%   - sig   Stress vector [n_el x 1]
%           sig(e) - Stress of bar e
%--------------------------------------------------------------------------
% It must provide as output:
%   - sig_cr   Critical stress vector [n_el x 1]
%              sig_cr(e) - Critical stress of bar e
%--------------------------------------------------------------------------

% Initialize
sig_cr = zeros(size(Tn,1),1);

% Iterate over all elements
for i = 1 : n_el
    % Find length vector components
    le_vec = x(Tn(i,2),:)-x(Tn(i,1),:);
    le = norm(le_vec);

    Inertia_e = pi/64*(mat(Tmat(i),5)^4-mat(Tmat(i),6)^4);
    
    sig_cr(i) = (pi^2*mat(Tmat(i),1)*Inertia_e)/(mat(Tmat(i),2)*le^2);

    if sig(i) < 0 && -sig(i) > sig_cr(i)
        fprintf('Bar %i will suffer buckling. \n', i);
    end

end

end