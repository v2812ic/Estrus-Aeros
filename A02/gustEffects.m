function [Lp,Dp,Tp,ax,az,Fi] = gustEffects(n,n_dof,n_el,n_i,m_p,mat,Tmat,x,cg,Tn,Aeronodes,Tnodes,g,L,D,T,Fw,Fext,gust)
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
%   - L     Previous Lift force
%   - D     Previous Drag force
%   - T     Previous Thrust force
%   - Fw    Weight vectors
%   - Fext  External forces vectors
%   - gust  Factor in which the gust increases aerodynamic loads
%--------------------------------------------------------------------------
% It must provide as output:
%   - Lp    New Lift force
%   - Dp    New Drag force
%   - Tp    New Thrust force
%   - ax    x-axis' induced acceleration
%   - az    z-axis' induced acceleration
%   - Fi    Inertial forces vectors
%--------------------------------------------------------------------------

if gust == 0
    % Say there is no gust
    fprintf('No gust here. Steady flight.\n');
    Lp = L; Dp = D; Tp = T; ax = 0; az = 0; Fi = zeros(n_dof,1);

else
    % Say there is gust
    fprintf('Gust of is present. Not in steady flight.\n');

    % Adjust aerodynamic loads
    Lp = gust*L;
    Dp = gust*D;
   
    % Declare the matrices
    A = zeros(3);
    B = zeros(3,1);
    C = zeros(3,1);
    
    % Find total weight
    W = sum(Fw+Fext);
    m = abs(W/g);

    % Define matrix A
    A(1,2) = m;
    A(2,1) = m;
    A(2,3) = -1;
    A(3,3) = 1/size(Tnodes,2)*sum(cg(Tnodes,3));

    % Define matrix B
    B(1) = Lp + W;
    B(2) = Dp;
    B(3) = Lp/size(Aeronodes,2)*sum(cg(Aeronodes,1)) - Dp/size(Aeronodes,2)*sum(cg(Aeronodes,3));

    % Solve the system
    C = A\B;
    ax = C(1);
    az = C(2);
    Tp = C(3);
    acc_vec = transpose([ax,0,az]);
    
    % Compute the mass of every node
    m_all_nodes = zeros(n,1);

    for i = 1 : n_el
        % Find the vector and length for each element
        le_vec = transpose(x(Tn(i,2),:)-x(Tn(i,1),:));
        le = norm(le_vec);
    
        % Calculate the mass of the bar
        m = mat(Tmat(i),3)*mat(Tmat(i),2)*le;
        m_all_nodes(Tn(i,1)) = m_all_nodes(Tn(i,1)) + m/2;
        m_all_nodes(Tn(i,2)) = m_all_nodes(Tn(i,2)) + m/2;
    
    end
    
    % Add passenger's mass
    m_all_nodes = m_all_nodes + transpose(m_p);

    % Initialize inertial forces vector
    Fi = zeros(n_dof,1);

    for i = 1 : n
        % Find the force for every node
        Fi(n_i*(i-1)+1:n_i*(i-1)+3) = -m_all_nodes(i)*acc_vec;
    end

end

end