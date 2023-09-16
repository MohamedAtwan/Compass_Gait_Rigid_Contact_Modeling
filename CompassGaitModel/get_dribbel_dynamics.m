clear;clc; close all;
syms ll lu mu ml Ju Jl q3 q5 q6 q3d q5d q6d g real
q = [q3 q5 q6];
qd = [q3d q5d q6d];
bi = sym('pi','real');
l(1) = Link([0 0 ll+lu 0]);
l(1).m = mu+ml;
l(1).r = [((ll+sym('0.35','real'))+ll/2)/2 0.0 0.0 ];
l(1).I = [0 0 0; 0 0 0; 0 0 (Ju+mu*((ll+sym('0.35','real'))-((ll+sym('0.35','real'))+ll/2)/2)^2+Jl+ml*(((ll+sym('0.35','real'))+ll/2)/2-ll/2)^2)/2];
l(1).G = 1;
l(1).Jm = 0;

l(2) = Link([0 0 lu 0]);
l(2).m = mu;
l(2).r = [lu-sym('0.35','real') 0.0 0.0 ];
l(2).I = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 Ju];
l(2).G = 1;
l(2).Jm = 0;

l(3) = Link([0 0 ll 0]);
l(3).m = ml;
l(3).r = [ll/2 0.0 0.0];
l(3).I = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 Jl];
l(3).G = 3;
l(3).Jm = 0;

r = SerialLink(l,'name','walker');
r.gravity = [-g 0 0];
% r.base = troty(bi/2)*trotz(-bi);
% J = r.jacob0([q3 -bi+q5 -q6]);
% R = roty(-pi/2);
% simplify(R*J(1:3,:))
tauG = r.gravload(q);
D = r.inertia(q);
C = r.coriolis(q,qd);
% r.links(1).dyn.r = [0.0 0.0 ((ll+sym('0.35'))+ll/2)/2];

% r.plot([pi/4 bi-bi/4 -pi/4]);