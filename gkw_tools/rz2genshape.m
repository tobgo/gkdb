% Computes the shape parametrisation to be used in the GKDB from R,Z FS description
%
   %  [alpha,beta,alpha_pr,beta_pr,R0,Z0]=rz2genshape(R,Z,r0,N,doplots);
%
% Inputs:
%  R,Z:   flux surfaces description in [m]
%         Size of the array assumed to be (nrho,npol), take care of it: no check performed!
%         The array is also assumed to have no double points
%  r0:    minor radius for which the parametrisation is performed
%  N:     Number of terms kept in the Fourier expansion
%  doplots: if 1 (default) perform plots to check the quality of the parametrisation
%           if -1, no plots at all
%           otherwise, no plots except if chi2 is to low  
%
% Outputs (as defined in IOGKDB.pdf):
%  alpha,beta:  cos and sin coefficients, in [m]
%  alpha_pr, beta_pr: radial derivative of alpha and beta
%
% YC - 01.10.2014

function [alpha,beta,alpha_pr,beta_pr,R0,Z0,err_out]=rz2genshape(R,Z,r0,N,doplots)

if ~exist('doplots')||isempty(doplots)
 doplots=1;
end 

if ~exist('N')||isempty(N)
 N=8;
end 

Nr=size(R,1);
Nth=size(R,2);
err_th = 0.003;

% minor radius
r = (max(R,[],2)-min(R,[],2))/2;

% reference point
R0dum = (max(R,[],2)+min(R,[],2))/2;
Z0dum = (max(Z,[],2)+min(Z,[],2))/2;

if r0<min(r) | r0>max(r)
 error('No extrapolation allowed, check value of r0')
end

R0 = interpos(2,r,R0dum,r0,-0.1);
Z0 = interpos(2,r,Z0dum,r0,-0.1);

% FS distance to the reference point
a = sqrt((R-R0).^2+(Z-Z0).^2);

% poloidal angle
th = atan2(-(Z-Z0)./a,(R-R0)./a);
th = mod(th+2*pi,2*pi);

% interpolate on a regular theta grid
th_grid = linspace(0,2*pi,Nth+1);
th_grid = th_grid(1:end-1);

[ath,coeffs]=deal(NaN.*zeros(Nr,Nth));

for ii=1:Nr
 [dum I]=sort(th(ii,:));
 ath(ii,:)=interpos(dum,a(ii,I),th_grid,-0.3);
 coeffs(ii,:)=fft(ath(ii,:));
end

alpha_all(:,1:N+1)  = real(coeffs(:,1:N+1))./Nth;
alpha_all(:,2:N+1)  = 2*alpha_all(:,2:N+1);
beta_all(:,1:N+1)   = -imag(coeffs(:,1:N+1))./Nth;
beta_all(:,2:N+1)   = 2*beta_all(:,2:N+1);


% check the parametrisation quality
C = cos([0:N]'*th_grid);
S = sin([0:N]'*th_grid);

a_out = repmat(alpha_all(:,1:N+1),[1 1 Nth]).*shiftdim(repmat(C,[1 1 Nr]),2) + ...
        repmat(beta_all(:,1:N+1),[1 1 Nth]).*shiftdim(repmat(S,[1 1 Nr]),2);
a_out = squeeze(sum(a_out,2));

R_out = R0 + a_out.*repmat(cos(th_grid),Nr,1);
Z_out = Z0 - a_out.*repmat(sin(th_grid),Nr,1);

err=sum(abs((ath-a_out)./ath),2)./Nth;

% Warning if the parametrisation quality is poor
dum = interpos(2,r,err,r0,0.);
if dum>=err_th;
  disp('Warning, parametrisation quality is poor at r=r0, consider increasing N')
end
err_out=dum;

if doplots==1

 figure
 plot(r./max(r),err)
 hold on
 plot([0 1],[err_th err_th],'r')
 xlabel('r/r_{max}')
 ylabel('relative error on the parametrisation')
 axis([0 1 0 0.05]) 

 figure
 plot(R',Z','b')
 hold on
 plot(R_out(err<err_th,:)',Z_out(err<err_th,:)','r--')

end

% interpolate the shaping coefficients at the desired location
for ii=1:N+1
 [alpha(ii) alpha_pr(ii)]=interpos(2,r,alpha_all(:,ii),r0,-0.1); 
 [beta(ii) beta_pr(ii)]=interpos(2,r,beta_all(:,ii),r0,-0.1); 
end



