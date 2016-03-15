clear
clear all
clc

syms phi theta p a b c phiStart phiEnd thetaStart thetaEnd

x = (p*sin(phi)*cos(theta))*a
y = (p*sin(phi)*sin(theta))*b
z = (p*cos(phi))*c



jacobian = [diff(x, p) diff(x, theta) diff(x, phi); diff(y, p) ...
    diff(y, theta) diff(y, phi); diff(z, p) diff(z, theta) diff(z, phi);]

% absolute value of the determinant (does not work with absolute value)
determ = -1 * det(jacobian)

% limits always from 0 to 1 for p
inner = int(determ, p, 0, 1);

% user specified limits phi
middle = int(inner, phi, phiStart, phiEnd);

% user specified limits theta
outer = int(middle, theta, thetaStart, thetaEnd)
