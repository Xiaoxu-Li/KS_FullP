function y = spherical_harmonic_xyz(l,m,x,y,z)
% function [real, imag] = spherical_harmonic_xyz(l,m,x,y,z)
% return Y_lm(x,y,z)

r = (x^2+y^2+z^2)^(1/2);
if r == 0
    error('error r=0');
else
    theta = acos(z/r);
    if y==0
        if x>=0
            phi = pi/2;
        else
            phi = -pi/2;
        end
    elseif y>0
        phi = atan(x/y);
    else
        if x>=0
            phi = atan(x/y)+pi;
        else
            phi = atan(x/y)-pi;
        end
    end
    [real, imag] = spherical_harmonic(l,m,theta,phi);
end
y = real + i*imag;
return;