% BEC with a harmonic oscillator potential of a stirrer corresponding a far-blue detuned Gaussian laser beam 2D

function v = V_gauss_2D(x,y)
r = (x^2+y^2)^(1/2);
v = 10*r^2+exp(-((x-1)^2+y^2));
return;
