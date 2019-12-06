% function y = V_hydrogen(r) 2D


function v = V_hydr_2D(x,y)
r = x^2+y^2;
v = -1.0/(r+5);
return; 