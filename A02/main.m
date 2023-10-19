%-------------------------------------------------------------------------%
% ASSIGNMENT 02
%-------------------------------------------------------------------------%
% Date: September 2023
% Author/s: Gerard Morales Riera, Víctor Rodríguez Romero
%
% Disclaimer: This code is prepared to solve a 3D structure. It might 
% not work for particular cases in 2D 
clear;
clc;
close all;

%% INPUT DATA [SI]

F = 800; %[N]
g = 9.81; %[N/kg]
D1 = 18e-3; %[m]
d1 = 7.5e-3; %[m]
D2 = 3e-3; %[m]
rho1 = 3350; %[kg/m^3]
rho2 = 950; %[kg/m^3]
Young1 = 75e9; %[Pa]
Young2 = 147e9; %[Pa]
thermal_coeff = 20e-6; %[K^-1]
Inertia = 1.2e-9; %[m^4]
Delta_T = 0; %[ºCelsius or Kelvin]
gust = 1.1; %[non-dim. factor]

%% PREPROCESS

% Nodal coordinates matrix creation
%  x(a,j) = coordinate of node a in the dimension j
x = [[0.85 -0.85/2 -0.9];
     [0.85 0.85/2 -0.9];
     [0.85 0 0];
     [-0.85 0 0];
     [-0.85 -3.2 0];
     [-0.85 3.2 0];
     [0 0 0]
    ];
 
% Connectivities matrix ceation
%  Tn(e,a) = global nodal number associated to node a of element e
Tn = [[1 2];
      [2 3];
      [1 3];
      [3 5];
      [3 6];
      [3 7];
      [5 7];
      [6 7];
      [4 5];
      [4 6];
      [4 7];
      [1 4];
      [2 4];
      [1 7];
      [2 7];
      [1 5];
      [2 6];
     ];

% External force matrix creation
%  Fdata(k,1) = node at which the force is applied
%  Fdata(k,2) = DOF (direction) at which the force is applied
%  Fdata(k,3) = force magnitude in the corresponding DOF
Fdata = [[1 3 -735.75];
         [2 3 -735.75];
        ];

% Additional mass for the passenger
m_p = [75,75,0,0,0,0,0];

% Fix nodes matrix creation
%  fixNod(k,1) = node at which some DOF is prescribed
%  fixNod(k,2) = DOF prescribed
%  fixNod(k,3) = prescribed displacement in the corresponding DOF (0 for fixed)
fixNod = [[7 3 0];
          [4 1 0];
          [4 2 0];
          [4 3 0];
          [5 1 0];
          [5 3 0];
          ];

% Aeronodes: nodes where aerodynamic forces (Lift, Drag) are applied
% Tnodes: nodes where thrust is applied
Aeronodes = [3, 4, 5, 6, 7];
Tnodes = [1, 2];

% Material data
%  mat(m,1) = Young modulus of material m
%  mat(m,2) = Section area of material m
%  --more columns can be added for additional material properties--
mat = [% Young M.   Section A.   thermal_coeff   Inertia
        [Young1, pi*(D1^2-d1^2)/4, rho1, thermal_coeff, D1, d1]; % M1
        [Young2,   pi*D2^2/4, rho2, thermal_coeff, D2, 0]; % M2
      ];

% Material connectivities
% Tmat(e) = Row in mat corresponding to the material associated to element e 
Tmat = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2];

%% SOLVER

% Dimensions
n_d = size(x,2);              % Number of dimensions 
n_i = n_d;                    % Number of DOFs for each node
n = size(x,1);                % Total number of nodes
n_dof = n_i*n;                % Total number of degrees of freedom
n_el = size(Tn,1);            % Total number of elements (bars)
n_nod = size(Tn,2);           % Number of nodes for each element (always 2)
n_el_dof = n_i*n_nod;         % Number of DOFs for each element 

% Computation of the DOFs connectivities
Td = connectDOFs(n_el,n_nod,n_i,n_el_dof,Tn);

% Computation of element stiffness matrices
Kel = computeKelBar(n_d,n_el,n_el_dof,x,Tn,mat,Tmat);

% Global matrix assembly
KG = assemblyKG(n_el,n_el_dof,n_dof,Td,Kel);

% Initial internal stress force
Fint_0 = initialStressForce(n_dof,n_el,n_i,x,Tn,mat,Tmat,Delta_T);

% Structure's weight computation
Fw = computeWeight(n_dof,n_el,mat,Tmat,x,Tn,g);

% Global force vector assembly
Fext = computeF(n_i,n_dof,Fdata);

% Find and establish reference from CG
cg = referenceFromCG(n,n_el,n_i,m_p,mat,Tmat,x,Tn);

% Lift and Drag computation 
[L,D,T] = computeLD(n,Fw,Fext,x,Aeronodes,Tnodes);

% Gust's effects on forces
[Lp,Dp,Tp,ax,az,Fi] = gustEffects(n,n_dof,n_el,n_i,m_p,mat,Tmat,x,cg,Tn,Aeronodes,Tnodes,g,L,D,T,Fw,Fext,gust);

% Total force computation
Ft = computeTF(Fext,Fint_0,Fw,Fi,Lp,Dp,Tp,Aeronodes,Tnodes);

% Apply conditions 
[vL,vR,uR] = applyCond(n_i,n_dof,fixNod);

% System resolution
[u,R] = solveSys(vL,vR,uR,KG,Ft);

% Compute strain and stresses
[eps,sig] = computeStrainStressBar(n_d,n_el,u,Td,x,Tn,mat,Tmat,Delta_T);

% Buckling analysis
[sig_cr] = analyseBuckling(n_el,x,Tn,mat,Tmat,sig);

%% POSTPROCESS

% Plot stress 3D
plotBarStress3D(x,Tn,u,sig,1);

%% PARTICULAR SCENARIO - B3

Linc = 0.1; %percentage of lift increase
Dinc = 0.1; %percentage of drag increase

% New Thrust Tp computation

% Acceleration computation
%% POSTPROCESS B3
plotBarStress3D(x,Tn,u,sig,1);