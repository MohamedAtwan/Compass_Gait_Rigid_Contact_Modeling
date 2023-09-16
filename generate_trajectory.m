[cx, cxd, cxdd, t] = cubic_traj(0, 5, 0, 0, 0.5176, 0);
[cy1, cy1d, cy1dd, t1] = cubic_traj(0, 2.5, 0, 0, 0.3, 0.1);
[cy2, cy2d, cy2dd, t2] = cubic_traj(2.5, 5, 0.3, 0.1, 0.0, 0.0);
cy = [cy1 cy2];
cyd = [cy1d cy2d];
cydd = [cy1dd cy2dd];

cxt = timeseries(cx, t);
cxdt = timeseries(cxd, t);
cxddt = timeseries(cxdd, t);

cyt = timeseries(cy, t);
cydt = timeseries(cyd, t);
cyddt = timeseries(cydd, t);



