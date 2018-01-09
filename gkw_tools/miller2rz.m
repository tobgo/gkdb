% Compute the R,Z description of a FS from the Miller parametrisation (GKW convention)
% 
%   [R,Z] = miller2rz(r,Rmil,Zmil,k,d,z,dRmildr,dZmildr,sk,sd,sz,Nth)
% 
% Inputs:
%  Miller parametrisation coeffficients
%
% Outputs:
%  R,Z: flux surface description in [m]
%       
%

function [R,Z] = miller2rz(r,Rmil,Zmil,k,d,z,dRmildr,dZmildr,sk,sd,sz,Nth)



if ~exist('Nth')||isempty(Nth)
 Nth=120;
end

th_grid = linspace(0,2*pi,Nth+1);
th_grid = th_grid(1:end-1);

Nr=3;
r0=r;

dr=[-0.01 ; 0 ; 0.01].*r0;
r_all = repmat(r + dr,1,Nth);
Rmil_all = repmat(Rmil + dRmildr.*dr,1,Nth);
Zmil_all = repmat(Zmil + dZmildr.*dr,1,Nth);
k_all = repmat(k + sk.*k./r.*dr,1,Nth);
d_all = repmat(d + sd.*sqrt(1-d.^2)./r.*dr,1,Nth);
z_all = repmat(z + sz./r.*dr,1,Nth);



R = Rmil_all + r_all.*cos(repmat(th_grid,Nr,1) + d_all.*sin(repmat(th_grid,Nr,1))); 
Z = Zmil_all + r_all.*k_all.*sin(repmat(th_grid,Nr,1)+z_all.*sin(2.*repmat(th_grid,Nr,1)));


