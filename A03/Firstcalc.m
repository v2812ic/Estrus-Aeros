function [A, Izz, C, lp] = Firstcalc(t1, t2, h1, h2, b, Me, g, eqmass1, mxlim1, eqmass2, mxlim2, eqlift1, lxlim1, eqlift2, lxlim2)
% FIRSTCALC 
% Computes the first geometrical data needed:
% Area of the cross section - A
% Cross inertia of the section - Izz
% Position of the centroid respect to the middle of the h1 side - cent
% Equilibrium parameter - lp

% It is being considered that the section is symmetric, therefore the y
% coordinate of the centroid is not necessary and some simpifications are
% being taken into account when calculating the area, such as neglecting
% the excess/defect of area in the joints.
% Also, the inertia is calculated considering slender beams.

ld = sqrt(b^2 + (h1-h2)^2/4); % longitude of the diagonal beam

A = t2*(h1+h2) + 2*t1*ld;
C = (t2*h2*b + b*t1*ld)/A;

beta = atan((h1-h2)/(2*b)); % angle of slope of the diagonal beam
yd = (h1+h2)/4; % gc y-coordinate of the diagonal beam [m]

% Cross inertia of the structure: lateral beams + diagonal beams + Steiner 
% on diagonal beams [m^4]
Izz = t2/12*(h1^3+h2^3) + t1/6*ld^3*(sin(beta))^2 + 2*t1*ld*yd^2; 

% Equilibrium parameter lp [N/m]
syms x l
Weight = (Me + int(eqmass1, x, mxlim1) + int(eqmass2, x, mxlim2))*g;
Lift = int(eqlift1, x, lxlim1) + int(eqlift2, x, lxlim2);
lp = double(solve(Weight == Lift, l));

end











