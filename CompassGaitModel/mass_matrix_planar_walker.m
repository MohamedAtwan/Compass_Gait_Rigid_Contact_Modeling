function M = mass_matrix_planar_walker(q,lu,ll,m1,m2,m3)

Jvc1 = [0 0 0; -(ll/2)*cos(q(1)) 0 0; -(ll/2)*sin(q(1)) 0 0];

Jvc2 = [0 0 0; (7/20)*cos(q(1)+q(2))-(ll+lu)*cos(q(1)) (7/20)*cos(q(1)+q(2)) 0;
        -(ll+lu)*sin(q(1))+(7/20)*sin(q(1)+q(2)) (7/20)*sin(q(1)+q(2)) 0];
Jvc3 = [0 0 0; -(lu+ll)*cos(q(3))+lu*cos(q(1)+q(2))+(ll/2)*cos(sum(q)) lu*cos(q(1)+q(2))+(ll/2)*cos(sum(q)) (ll/2)*cos(sum(q));
        -(lu+ll)*sin(q(1))+lu*sin(q(1)+q(2))+(ll/2)*sin(sum(q)) lu*sin(q(1)+q(2))+(ll/2)*sin(sum(q)) (ll/2)*sin(sum(q))];
D = m1*(Jvc1'*Jvc1)+m2*(Jvc2'*Jvc2)+m3*(Jvc3'*Jvc3);
M=D;