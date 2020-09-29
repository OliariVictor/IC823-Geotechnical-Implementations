% Viga em fundacao elastica
% Hipotese de Winckler
clear all
close all
clc
Lx = 2.0; nx = 20;
Ly = 1.0; ny = 10;
e_soil = 1; poisson = 0.4; p = 1.0;
n = nx*ny;
L = Lx/nx;
B = Ly/ny;

x0 = L/2:L:Lx-L/2;
y0 = B/2:B:Ly-B/2;

[xc, yc] = meshgrid(x0,y0);
xc = reshape(xc, [n,1]);
yc = reshape(yc, [n,1]);

infl = zeros(n,n);
for i=1:n
    for j=1:n
        xa = abs(xc(i) - xc(j)) - L/2;
        ya = abs(yc(i) - yc(j)) - B/2;
        rx = sqrt((L + xa)^2 + ya^2);
        ry = sqrt((B + ya)^2 + xa^2);
        ra = sqrt((xa^2 + ya^2));
        rc = sqrt((L + xa)^2 + (B + ya)^2);
        infl(i,j) = 1/pi * (xa*log((ya+ra)*(L+ya+rc)/(ya+rx)/(B+ya+ry))+...
            ya*log((xa+ra)*(xa+ry))+L*log((B+ya+rc)/(ya+rx))+...
            B*log((L+xa+rc)/(xa+ry))+L*log((B+ya+rc)/(ya+rx)));
    end
end
bv = [zeros(n,1);p];
k = [infl,-ones(n,1);L*B*ones(1,n),0];
d=k\bv;
t=d(1:n);
dr=d(n+1);
w=dr*(1-poisson^2)/e_soil
