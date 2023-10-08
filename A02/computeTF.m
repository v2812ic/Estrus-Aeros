function Ft = computeTF(Fext, Fint_0, Fw, L, D, T, Aeronodes, Tnodes)
%--------------------------------------------------------------------------
% The function takes as inputs:
%--------------------------------------------------------------------------
% It must provide as output:
%--------------------------------------------------------------------------

Ft = Fext + Fint_0 + Fw;
Ft(Aeronodes*3) = Ft(Aeronodes*3) + L/size(Aeronodes, 2);
Ft(Aeronodes*3 - 2) = Ft(Aeronodes*3 - 2) + D/size(Aeronodes, 2);
Ft(Tnodes*3 - 2) = Ft(Tnodes*3 - 2) + T/size(Tnodes, 2);
end
