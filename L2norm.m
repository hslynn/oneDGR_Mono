function l2 = L2norm(u)
Globals1D;
l2 = sqrt(sum(u(1:end-1).*u(1:end-1)*dx));
return
