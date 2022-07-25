function [pPlus, pMinus] = deri(u)
Globals1D;
pPlus = zeros(1, size(u)(2));
pMinus = zeros(1, size(u)(2));
pPlus(1:end-1) = (u(2:end) - u(1:end-1))./dx;
pMinus(2:end) = (u(2:end) - u(1:end-1))./dx;

pPlus(end) = (u(end) - u(end))./dx;
pMinus(1) = (u(1) - u(1))./dx;
return

