function xc = mostSolicitatedPosition(t1,t2,b,L1,L2,h1,h2,C,Nel,Izz,Fy,Mz)

% Initial values
l_div = (L1+L2)/Nel; % length of each division
ymax = max(h1,h2)/2; % max height
Ain = b*(h1+h2)/2; % internal area of the cross section
m = (h1/2-h2/2)/b; % slope of upper and lower sides
sig_eq_max = 0;
syms s;

% Constant integrals
    % Upper side
    y1(s) = h2/2 + s*(h1/2-h2/2)/b;
    ty1(s) = t1*y1;
    g1(s) = int(ty1,s,0,s);
    r1 = m^2/(1+m^2)*(h1/2-m*C);

    % Left side
    y2(s) = h1/2 - s;
    ty2(s) = t2*y2;
    g2(s) = int(ty2,s,0,s);
    r2 = C;

    % Lower side
    y3(s) = -h1/2 + s*(h1/2-h2/2)/b;
    ty3(s) = t1*y3;
    g3(s) = int(ty3,s,0,s);
    r3 = r1;

    % Right side
    y4(s) = -h2/2 + s;
    ty4(s) = t2*y4;
    g4(s) = int(ty4,s,0,s);
    r4 = b - C;

for i = 1 : Nel
    i
    % Find the normal stress due to Mz
    sig_max = -(Mz(i,1)*ymax)/Izz

    % Find the shear stress due to Fy
    qo_1(s) = -(Fy(i,1))/Izz*g1(s);
    qo_2(s) = -(Fy(i,1))/Izz*g2(s) + qo_1(sqrt(b^2+(h1/2-h2/2)^2));
    qo_3(s) = -(Fy(i,1))/Izz*g3(s) + qo_2(h1);
    qo_4(s) = -(Fy(i,1))/Izz*g4(s) + qo_3(sqrt(b^2+(h1/2-h2/2)^2));

    q0 = -1/(Ain*2)*(r1*int(qo_1,s,0,sqrt(b^2+(h1/2-h2/2)^2))+r2*int(qo_2,s,0,h1)+r3*int(qo_3,s,0,sqrt(b^2+(h1/2-h2/2)^2))+r4*int(qo_4,s,0,h2));

    q1(s) = qo_1(s) + q0;
    q2(s) = qo_2(s) + q0;
    q3(s) = qo_3(s) + q0;
    q4(s) = qo_4(s) + q0;

    tau1(s) = q1(s)/t1;
    tau2(s) = q2(s)/t2;
    tau3(s) = q3(s)/t1;
    tau4(s) = q4(s)/t2;

    tau_max = tau3(0);

    sig_eq = sqrt(sig_max^2+3*tau_max^2);

    if sig_eq > sig_eq_max
        nc = i;
        sig_eq_max = sig_eq;

    end
    
    previous = sig_eq;

end

xc = (nc-1)*l_div;
fprintf('Most critical point at %2.2f m with sig: %4.4f Pa.\n',xc,sig_eq_max);

end