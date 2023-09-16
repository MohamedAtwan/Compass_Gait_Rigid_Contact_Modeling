% close all;
lu = 0.47;
ll = 0.43;
mu = 3; ml = 1.6;
Ju = 0.108; Jl = 0.059;
g = 9.81;
l(1) = Link([0 0 ll+lu 0]);
l(1).m = mu+ml;
l(1).r = [((ll+0.35)+ll/2)/2 0.0 0.0 ];
l(1).I = [0 0 0; 0 0 0; 0 0 (Ju+mu*((ll+0.35)-((ll+0.35)+ll/2)/2)^2+Jl+ml*(((ll+0.35)+ll/2)/2-ll/2)^2)/2];
l(1).G = 1;
l(1).Jm = 0;

l(2) = Link([0 0 lu 0]);
l(2).m = mu;
l(2).r = [lu-0.35 0.0 0.0 ];
l(2).I = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 Ju];
l(2).G = 1;
l(2).Jm = 0;

l(3) = Link([0 0 ll 0]);
l(3).m = ml;
l(3).r = [ll/2 0.0 0.0];
l(3).I = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 Jl];
l(3).G = 1;
l(3).Jm = 0;

r = SerialLink(l,'name','walker');
r.gravity = [-g 0 0];
% r.base = troty(pi/2)*trotz(-pi);
r.plot([0 -pi+pi/4 -0]);