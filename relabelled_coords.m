function [q_new, jacc] = relabelled_coords(q)
q_new  = zeros(size(q));
q_new(1) = q(1)+q(2);
q_new(2) = -q(2);

jacc = [1 0; 1 -1];
