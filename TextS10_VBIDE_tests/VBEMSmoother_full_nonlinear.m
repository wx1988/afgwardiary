function [Estinfo] = VBEMSmoother_full_nonlinear(spikes,Constants,Initinfo,FieldMatrices,Basis,Thetainfo)
% -------------------------------------------------------------
% Function VBEMFilter_full_nonlinear
% inputs:
%   spikes, list of structs of coordinate
%   Constants, some fixed paramter, quite a lot of data here
%   InitInfo, I think this filed is mainly related with the so-called state
%   FieldMatrices, seems to be the state-space model relation. 
%   Basis, the basis function
%   Thetainfo, the parameter part. 
%
% outputs:
%   Estinfo, list of estimated data, each item for one week.
%
% -------------------------------------------------------------


PSIxinv = inv(FieldMatrices.PSIx);
SigmaW = inv(Thetainfo.Meanprecmat);
PSIxinvV = PSIxinv*FieldMatrices.V(:,:,1);

rho = Thetainfo.rho;
varrho = Thetainfo.varrho;
Aest = rho*Thetainfo.kernelpar*PSIxinvV;
beta = Constants.beta;

xestprior = zeros(Basis.nx,Constants.N);
xestpost = zeros(Basis.nx,Constants.N);
xbeta = zeros(Basis.nx,Constants.N);
xestRTS = zeros(Basis.nx,Constants.N);

Sigmaprior = zeros(Basis.nx,Basis.nx,Constants.N);
Sigmapost = zeros(Basis.nx,Basis.nx,Constants.N);
Sigmabeta = zeros(Basis.nx,Basis.nx,Constants.N);
SigmaRTS = zeros(Basis.nx,Basis.nx,Constants.N);

% TODO, why fixed constant 50?
Sigmaprior(:,:,1) = 50*eye(Basis.nx);
xestprior(:,1) = Initinfo(1).xestprior;
Sigmapost(:,:,1) = 30*eye(Basis.nx);
xestpost(:,1) = Initinfo(1).xestpost;

mu = Thetainfo.muest;
theta = Thetainfo.thetaest;
Qinv = inv(SigmaW);
AQinvA = (rho^2 + varrho)*(Thetainfo.kernelpar^2 + Thetainfo.varkernelpar)*...
    PSIxinvV'*Qinv*PSIxinvV;
dt = Constants.dt;
y = reshape(spikes,[],size(spikes,3));
Estinfo = Initinfo;

options = foptions;
options(14) = 2000;
% options(9) = 1;
options(2) = 0.1;
options(3) = 0.1;

for i = 2:Constants.N
    xestprior(:,i) = Aest*xestpost(:,i-1) + Constants.dt*theta;
    Sigmaprior(:,:,i) = Aest*Sigmapost(:,:,i-1)*Aest' + SigmaW;
    
    Sigmatilde = inv(inv(Sigmapost(:,:,i-1)) + AQinvA);
    Sigmastar = inv(Qinv - Qinv*Aest*Sigmatilde*Aest'*Qinv);
    mustar = Sigmastar*(Qinv*Aest*Sigmatilde*(inv(Sigmapost(:,:,i-1))*xestpost(:,i-1) - Aest'*Qinv*theta*dt) + Qinv*theta*dt);
    
    %Scan data for NaNs
    spikecoords = spikes(i).Coords;
    
    if ~isempty(spikecoords)
        phieval = zeros(Basis.nx,size(spikecoords,1));
        for j = 1:Basis.nx
            phieval(j,:) = LocalisedKernelPhi_Cont(spikecoords(:,1),spikecoords(:,2),Basis.mu1(j),Basis.mu2(j),Basis.tau1(j),Basis.tau2(j))';
        end
    else phieval = zeros(Basis.nx,size(spikecoords,1));
    end
    
    myint = @(xx)  SingleIntForward(xx,Constants,Basis,Thetainfo);
    myint2 = @(xx)  MultiIntForward(xx,Constants,Basis,Thetainfo);
    f = @(xx) -(sum(mu + beta*phieval'*xx') - myint(xx') - ((xx' - mustar)'/Sigmastar)*(xx' - mustar)./2);
    gradf = @(xx) -(sum(beta*phieval,2)' - myint2(xx')' - xx/Sigmastar + mustar'/Sigmastar);
    %[temp,options,t1,t2,t3] = scg(f,xestprior(:,i)',options,gradf);
    [xestpost(:,i)] = scg(f,xestprior(:,i)',options,gradf)';
    myint3 = @(xx) MultiIntForward2(xx,Constants,Basis,Thetainfo);
    temp = inv(Sigmastar) + myint3(xestpost(:,i));
    Sigmapost(:,:,i) = inv(temp);
    [i max(xestpost(:,i)) min(xestpost(:,i))]
    Estinfo(i).xestpost = xestpost(:,i);
    Estinfo(i).xestprior = xestprior(:,i);
    Estinfo(i).PKalman = Sigmapost(:,:,i);
end


Sigmabeta(:,:,i) = 9*eye(Basis.nx);
xbeta(:,end) = xestpost(:,end);

for i = Constants.N-1:-1:1
    
    spikecoords = spikes(i+1).Coords;
    if ~isempty(spikecoords)
        phieval = zeros(Basis.nx,size(spikecoords,1));
        for j = 1:Basis.nx
            phieval(j,:) = LocalisedKernelPhi_Cont(spikecoords(:,1),spikecoords(:,2),Basis.mu1(j),Basis.mu2(j),Basis.tau1(j),Basis.tau2(j))';
        end
    else phieval = zeros(Basis.nx,size(spikecoords,1));
    end
    
    %     Hacky way linearizing around filtered estimated
    %     Sigmadash = inv(inv(Sigmabeta(:,:,i+1)) + myint3(xestpost(:,i+1)));
    %     mudash = xestpost(:,i+1) + Sigmadash*(inv(Sigmabeta(:,:,i+1))*(xbeta(:,i+1) - xestpost(:,i+1)) + sum(beta*phieval,2) - myint2(xestpost(:,i+1)));
    
    %   Proper way
    f = @(xx) -(sum(mu + beta*phieval'*xx') - myint(xx') - ((xx' - xbeta(:,i+1))'/Sigmabeta(:,:,i+1))*(xx' - xbeta(:,i+1))./2);
    gradf = @(xx) -(sum(beta*phieval,2)' - myint2(xx')' - xx/Sigmabeta(:,:,i+1) + xbeta(:,i+1)'/Sigmabeta(:,:,i+1));
    mudash = scg(f,xestpost(:,i+1)',options,gradf)';
    Sigmadash = inv(inv(Sigmabeta(:,:,i+1)) + myint3(mudash));
    
    Sigmatilde = inv(inv(Sigmadash) + Qinv);
    Sigmabeta(:,:,i) = inv(AQinvA - Aest'*Qinv*Sigmatilde*Qinv*Aest);
    xbeta(:,i) = Sigmabeta(:,:,i)*(-dt*Aest'*Qinv*theta + Aest'*Qinv*Sigmatilde*(inv(Sigmadash)*mudash + Qinv*dt*theta));
    
    SigmaRTS(:,:,i) = inv(inv(Sigmapost(:,:,i)) + inv(Sigmabeta(:,:,i)));
    xestRTS(:,i) = SigmaRTS(:,:,i)*(inv(Sigmapost(:,:,i))*xestpost(:,i) + inv(Sigmabeta(:,:,i))*xbeta(:,i));
    
    i
    
    Estinfo(i).xestRTS = xestRTS(:,i);
    Estinfo(i).PRTS = SigmaRTS(:,:,i);
    
    temp = Qinv + inv(Sigmabeta(:,:,i+1))+ myint3(xestRTS(:,i+1));
    Sigmatilde = inv(inv(Sigmapost(:,:,i)) + Qinv);
    Crosscov = Sigmatilde*(Qinv)*inv(temp -  Qinv*Sigmatilde*Qinv);
    
    Estinfo(i).W = xestRTS(:,i)*xestRTS(:,i)' + SigmaRTS(:,:,i);
    Estinfo(i).S = xestRTS(:,i)*xestRTS(:,i+1)' + Crosscov;
    
end
%Update posterior for next time step
xestpost(:,1) = Aest*xestRTS(:,2) - theta*dt; %initial state to the next iteration
% Sigmapost(:,:,1) = SigmaW*inv(eye(Basis.nx)-A*A'); %initial state

