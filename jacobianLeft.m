function jacc = jacobianLeft(q)
jacc = [cos(sum(q))-cos(q(1)) cos(sum(q));
        sin(sum(q))-sin(q(1)) sin(sum(q))];