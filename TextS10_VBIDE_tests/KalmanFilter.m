function [Estinfo] = KalmanFilter(y,Initinfo,Constants,FieldMatrices,Basis)
% -------------------------------------------------------------
% Function KalmanFilter
% inputs:
% outputs:
% Description: In conjunction with Diggle estimator to put into account the
% system dynamics
%
% http://www.cs.unc.edu/~welch/kalman/index.html#Anchor-Rudolph-6296
% 
% -------------------------------------------------------------


SigmaW = 0.2^2*FieldMatrices.W;
SigmaR = FieldMatrices.R;

xestpost = zeros(Basis.nx,Constants.N);
sigma2estpost = zeros(Basis.nx,Basis.nx,Constants.N);
xestprior = zeros(Basis.nx,Constants.N);
sigma2estprior = zeros(Basis.nx,Basis.nx,Constants.N);

sigma2estprior(:,:,1) = 50*eye(Basis.nx);
xestprior(:,1) = Initinfo(1).xestprior;
xestpost(:,1) = Initinfo(1).xestpost;

Estinfo = Initinfo;
C = eye(Basis.nx);

S = C*sigma2estprior(:,:,1)*C' + SigmaR;
K = sigma2estprior(:,:,1)*C'*inv(S);
xestpost(:,1) = xestprior(:,1) + K*(y(:,1) - C*xestprior(:,1));
sigma2estpost(:,:,1) = (eye(Basis.nx) - K*C)*sigma2estprior(:,:,1);
Estinfo(1).xestpost = xestpost(:,1);