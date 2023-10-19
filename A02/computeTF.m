function Ft = computeTF(Fext,Fint_0,Fw,Fi,Lp,Dp,Tp,Aeronodes,Tnodes)
%--------------------------------------------------------------------------
% The function takes as inputs:
%--------------------------------------------------------------------------
% It must provide as output:
%--------------------------------------------------------------------------

Ft = Fext + Fint_0 + Fw + Fi;
Ft(Aeronodes*3) = Ft(Aeronodes*3) + Lp/size(Aeronodes,2);
Ft(Aeronodes*3-2) = Ft(Aeronodes*3-2) + Dp/size(Aeronodes,2);
Ft(Tnodes*3-2) = Ft(Tnodes*3-2) + Tp/size(Tnodes,2);

end
