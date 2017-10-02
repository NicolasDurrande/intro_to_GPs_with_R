function [uxi,uyi,uzi]=mogi_3D(G,nu,xs,ys,zs,a,p,xi,yi,zi)

% MOGI(G,nu,xs,ys,zs,a,p,xi,yi,zi) compute surface displacements and tilts created by
% a point source located beneath a topography. To account for toography, a
% first order solution in which the actual source to ground surface point
% is taken into account
% 
% [uxi,uyi,uzi]=mogi(G,nu,xs,ys,zs,a,p,xi,yi,zi)
% Parameters are 
% G = shear modulus in MPa, G = E/2(1+nu)
% nu = Poisson's ratio
% xs, ys, zs = source position (z axis is positive upward),
% a = source radius, p = source overpressure in MPa, 
% xi, yi, zi = location of ground surface points
%
% V. Cayol, LMV, sept 2017

DV=pi*a^3*p/G;
C=(1-nu)*DV/pi;

r = sqrt((xi-xs).^2+(yi-ys).^2);
f = r.^2+(zi-zs).^2;
uzi = C*(zi-zs)./(f.^(3/2));
ur = C*r./(f.^(3/2));
theta = atan2(yi-ys,xi-xs);
uxi = ur.*cos(theta);
uyi = ur.*sin(theta);

figure
subplot(1,3,1), pcolor(xi/1000,yi/1000,uxi), shading interp, xlabel('x(km)'),ylabel('y(km)'), colorbar
subplot(1,3,2), pcolor(xi/1000,yi/1000,uyi), shading interp, xlabel('x(km)'),ylabel('y(km)'), colorbar
subplot(1,3,3), pcolor(xi/1000,yi/1000,uzi), shading interp, xlabel('x(km)'),ylabel('y(km)'), colorbar