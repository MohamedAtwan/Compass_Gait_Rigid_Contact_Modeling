function Xdd = unconstrained_acceleration(q,qd,qdd)

A = [cos(sum(q))-cos(q(1)) cos(sum(q));
        sin(sum(q))-sin(q(1)) sin(sum(q))];
A_dot = [-sum(qd)*sin(sum(q))+qd(1)*sin(q(1)) -sin(sum(q))*sum(qd);
         sum(qd)*cos(sum(q))-qd(1)*cos(q(1)) cos(sum(q))*sum(qd)];
Xdd = A*qdd+ A_dot*qd;