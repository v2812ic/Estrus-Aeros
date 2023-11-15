%-------------------------------------------------------------------------%
% ASSIGNMENT 03 - (A)
%-------------------------------------------------------------------------%
% Date: November 2023
% Author/s: Morales G. // Rodríguez V.

clear;
close all;

%% INPUT DATA

% Material properties
E = 85e9; % 85GPa
% Cross-section parameters [m]
t1 = 1.5e-3;
t2 = 4e-3;
h1 = 0.5;
h2 = 0.25;
b = 0.775;
% Other data [IS]
g = 9.81;
L1 = 5;
L2 = 10;
L = L1+L2;
Me = 2550;
M = 35e3;

% Distribution of mass and lift plus limits
syms x l 

eqmass1 = M/(4*L) + 3*M/(2*L2^2)*(L1-x);
mxlim1 = [0, L1];

eqmass2 = M/(4*L);
mxlim2 = [L1, L];

eqlift1 = l*(0.8 -0.2*cos(pi*x/L1));
lxlim1 = [0, L1];

eqlift2 = l*(1-(x-L1)/L2)*(1+(x-L1)/L2);
lxlim2 = [L1, L];

% Number of elements for each part
nel = [3, 6, 12, 24, 48, 96];

%% PRECOMPUTATIONS

%First calculations: 
[A, Izz, C, lp] = Firstcalc(t1, t2, h1, h2, b, Me, g, eqmass1, mxlim1, eqmass2, mxlim2, eqlift1, lxlim1, eqlift2, lxlim2);

% Plot analytical solution
% fig = plotBeamsInitialize(L);

% Loop through each of the divisions
for k = 1:length(nel)

    %% PREPROCESS
    
    % Nodal coordinates
    % x(a,j) = coordinate of node a in the dimension j (there's only one
    % dimension)
    % n_k: number of beams in the subdivision
    n_k = nel(k);
    x = linspace(0, L, n_k+1);
    
    % Nodal connectivities // This vectors differs from the original alg.
    %  Tnod(e,a) = Coordinate (x) associated to node a of element e
    Tnod = zeros(n_k, 2);
    Tnod(:, 1) = x(1, 1:(length(x)-1));
    Tnod(:, 2) = x(1, 2:length(x));
    
    % Material properties matrix
    %  mat(m,1) = Young modulus of material m
    %  mat(m,2) = Section area of material m
    %  mat(m,3) = Section inertia of material m
    mat = [
        E, A, Izz;  % Material (1)
    ];

    % Material connectivities
    %  Tmat(e) = Row in mat corresponding to the material associated to element e 
    Tmat = ones(1, n_k);
       
    %% SOLVER
    
    % Computation of element stiffness matrices
    Kel = computeKelBar(n_k, Tnod, mat, Tmat);
    
    % Computation of the element force vector
    Fel = computeForce(Tnod, eqmass1, eqmass2, eqlift1, eqlift2, lp, L1);
   
    % Global Stiffness matrix
    
    % Force vector assembly
    
    % Compute:
    % u  - Displacements and rotations vector [ndof x 1]
    % pu - Polynomial coefficients for displacements for each element [nel x 4]
    % pt - Polynomial coefficients for rotations for each element [nel x 3]
    % Fy - Internal shear force at each elements's nodes [nel x nne]
    % Mz - Internal bending moment at each elements's nodes [nel x nne]
    
    %% POSTPROCESS
    
    % Number of subdivisions and plots
    nsub = Nel/nel(k);
    plotBeams1D(fig,x,Tnod,nsub,pu,pt,Fy,Mz)
    drawnow;
    
end

% Add figure legends
figure(fig)
legend(strcat('N=',cellstr(string(nel))),'location','northeast');