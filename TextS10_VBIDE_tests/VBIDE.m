%--------------------------------------------------------------------------
% Program VBEM for IDE
% Authors: Andrew Zammit Mangion
% Date: March 19 2012
% 
% Details: An offline VBEM algorithm for state (field) + parameter
% inference of a spatiotemporal state-space system with point process observations
%
% Input: filename (output filename)
%        Estprec ('y' or 'n'. If 'n' the identity matrix is used for the precision)
%--------------------------------------------------------------------------

function y = VBIDE(filename,Estprec)

load('Test_data');  % The Test_data might be a simulation data. 


Esttheta = 'y';
Estkernel = 'n';
Estrho = 'n';
Estmu = 'n';

%Refine space
J = 101;
s = linspace(s(1),s(end),J); %#ok<NODEF>
ds = (max(s)-min(s))/(J-1);
% [tau,cutoff]=FrequencyAnalysis(r,gest,s);         %length of basis functions 0.15 for Fine data
tau = 1.4175;   % TODO, how is this tau determined? The tau in S6 is 1.8, why diff?
cutoff= 0.18;   % TODO, what is this cutoff for? the frequency?
smax = s(end);
smin = 0;
[s1,s2] = meshgrid(s,s);

%Arrange in struct for par. passing
Constants.dt = dt;  % load from test_data.mat, 1 
Constants.ds = ds;  % the sample width in spatial
Constants.N = N;    % might be number of weeks. 
Constants.t = t;    % week index
Constants.s = s;    % 1*25, ds*(i-1)
Constants.s1 = s1;  % horizontal index
Constants.s2 = s2;  % vertical index
Constants.smax = smax;  % max value in spatial
Constants.smin = smin;  % min value in spatial
Constants.J = J;    % the sample number


%-------------------------------------------------------------------
% Field setup
%-------------------------------------------------------------------

%Field noise
sigmaW = 0.2;
sigmaR = 2;

spacing = 1/(2*cutoff*1.2);
mu1 = linspace(s(1),s(end),(s(end)-s(1))/spacing + 1);
mu2 = linspace(s(1),s(end),(s(end)-s(1))/spacing + 1);
Basis.nx = length(mu1)^2; % TODO, assuming the same discretization for the two axes

% Field basis functions
[C1] = meshgrid(mu1,mu2);
Basis.mu1 = reshape(C1,1,[]);
Basis.mu2 = reshape(C1',1,[]);
tau1(1:Basis.nx) = tau;
tau2(1:Basis.nx) = tau;
Basis.tau1 = tau1;
Basis.tau2 = tau2;
Basis.phi = LocalisedKernelPhi(s1,s2,Basis.mu1,Basis.mu2,tau1,tau2);
close all

for i =1:Basis.nx
    surf(Basis.phi(:,:,i)); shading interp; hold on
end

if Estkernel == 'y'
    % TODO, this is set to 'n' for now.
    % Kernel basis functions
    gaussf = @(s1,s2,mu1,mu2,sigma21,sigma22) exp(-(s1-mu1).^2/sigma21-(s2-mu2).^2/sigma22);
    KernelBasis.d = 1;
    KernelBasis.phi(:,:,1) = gaussf(s1,s2,IDEKernel.mu1,IDEKernel.mu2,IDEKernel.sigma2,IDEKernel.sigma2);
end

% Required matrices
% TODO, not clear the inner operation
Basis.inner = gaussinner(s1,s2,Basis.phi);

% This is to change the matrix into vector, this migh make it easy for
% calculation
Basis.Basisvec = zeros(J^2,Basis.nx);
for i = 1:Basis.nx
    Basis.Basisvec(:,i) = reshape(Basis.phi(:,:,i),[],1);
end


if Estkernel == 'y'
    % TODO, Estkernel is set to 'n' now.
    KernelBasis.PHI = zeros(length(s),length(s),KernelBasis.d,Basis.nx);
    for j = 1:KernelBasis.d
        for i = 1:Basis.nx
            KernelBasis.PHI(:,:,j,i) = conv2(Basis.phi(:,:,i),KernelBasis.phi(:,:,j),'same')*ds^2;
        end
    end
end


%Noise mean
b(1,1:Basis.nx) = 0;
Constants.d = length(b);    % just Basic.nx, 

%Setting up matrices and state-space model.
FieldMatrices.PSIx = Basis.inner; 
% TODO, this is 64*64 matrix, similar with covariance. 
% TODO, how is this related to inferenc?


% Function definitions
% TODO, Growth is 25*25 matrix, seems to be the raw data to generated
% simulated data
[S1x,S2x] = meshgrid(linspace(s(1),s(end),length(Growth)));

gaussf = @(s1,s2,mu1,mu2,sigma21,sigma22) exp(-(s1-mu1).^2/sigma21-(s2-mu2).^2/sigma22);
Growth_interp = @(s1,s2) interp2(S1x,S2x,Growth,s1,s2);
Avariance_cont =  @(s1,s2) interp2(S1x,S2x,Avariance_map,s1,s2); %#ok<NODEF>

%Exact reduced noise covariance matrix
NoiseKernel.K = gaussf(s1,s2,NoiseKernel.mu1,NoiseKernel.mu2,lnoise,lnoise); %#ok<NODEF>
Avariance_map = Avariance_cont(s1,s2);
Qphi = zeros(J,J,Basis.nx);
for i = 1:Basis.nx
    Qphi(:,:,i) = conv2(Avariance_map.*Basis.phi(:,:,i),NoiseKernel.K,'same')*ds^2;
end
FieldMatrices.Qn = gaussinner2(s1,s2,repmat(Avariance_map,[1,1,Basis.nx]).*Basis.phi,Qphi);
W = inv(FieldMatrices.PSIx)*FieldMatrices.Qn*inv(FieldMatrices.PSIx);
W2 = triu(W) + triu(W)' - diag(diag(W)); %Ensure symmetry!! (numerical errors)
FieldMatrices.W = W2;

if Estkernel == 'y' 
    FieldMatrices.V = FindV(Constants.s1,Constants.s2,Basis,KernelBasis);
end

%Growth vector
% TODO, seems to generated the ground truth here
Basis_vec_mat = reshape(Basis.phi,[],Basis.nx);
thetatrue = inv(Basis_vec_mat'*Basis_vec_mat)*Basis_vec_mat'*reshape(Growth_interp(s1,s2),[],1);


% --------------------
% Initialise variables
% --------------------
numiters = 50;

% -------------------------------------------------------
% Initialise field estimate
% -------------------------------------------------------
Estinfo(1).xestpost = zeros(Basis.nx,1);
Estinfo(1).xestprior = zeros(Basis.nx,1);
Estinfo(1).xestRTS = zeros(Basis.nx,1);
Estinfo(1).PKalman = zeros(Basis.nx,Basis.nx);
Estinfo(1).PRTS = zeros(Basis.nx,Basis.nx);
Estinfo = repmat(Estinfo(1),1,N);

% -------------------------------------------------------
% Initialise parameter estimate
% -------------------------------------------------------

% TODO, it seems that the initial variance on theta could be estiamted.
% TODO, why is this theta equals to the number of basis function?
Thetainfo(1).thetaest = zeros(Constants.d,1);
% Thetainfo(1).thetaest = thetatrue;
Thetainfo(1).vartheta = 1*eye(Constants.d);%1000*eye(Constants.d);

% TODO, why only one dimention for mu, 
% is this mu the density of all position before the transformation?
Thetainfo(1).muest = 1;
Thetainfo(1).varmu = 0.1;

% Is this precision generated based on the variance of the noise?
Thetainfo(1).precisionest = repmat(1/(sigmaW^2),Basis.nx,1);
Thetainfo(1).varprecision = repmat(1,Basis.nx,1);

Thetainfo(1).Meanprecmat = eye(Basis.nx);
% Thetainfo(1).Meanprecmat = inv(FieldMatrices.W);
Thetainfo(1).kernelpar = 0.001;
Thetainfo(1).varkernelpar = 0.1;

Thetainfo(1).rho = 0.5;
Thetainfo(1).varrho = 0.5;

% TODO, it seems that if the parameter is to be estimated, 
% There will be a variance associated, which might be for the iterative
% calculation

if Estmu == 'n'
    Thetainfo(1).muest = bias; Thetainfo(1).varmu = 0;
end
if Estrho == 'n'
    Thetainfo(1).rho = rho; Thetainfo(1).varrho = 0;
end
if Esttheta == 'n'
    Thetainfo(1).thetaest = thetatrue; Thetainfo(1).vartheta = 0*eye(Constants.d);
end
if Estkernel == 'n'
    Thetainfo(1).kernelpar = 1;
    Thetainfo(1).varkernelpar = 0;
    FieldMatrices.V(:,:,1) = FieldMatrices.PSIx;
end
if Estprec == 'n'
%     my_sum(i) = 0; for i = 1:N my_sum(i) = length(spikeaccept(i).Coords(:,1)); end
%     Educated_var_guess = var(log(my_sum/s(end).^2));
    Thetainfo(1).Meanprecmat = eye(Basis.nx);
end

Thetainfo = repmat(Thetainfo(1),1,numiters);
Constants.theta = b;

% -------------------------------------------------------
% Point process parameters
% -------------------------------------------------------
beta = 1;
Constants.beta = beta;


%-----------------------------------------------
% Inference
%-----------------------------------------------
%Estimate initial condition
EstIntensity = DiggleContSpace(spikeaccept(1).Coords,Constants);
y = Regress(EstIntensity,Basis,Thetainfo(1).muest ,beta);
FieldMatrices.R = sigmaR^2*eye(Basis.nx);
Estinfo(1) = KalmanFilter(y,Estinfo(1),Constants,FieldMatrices,Basis);

% Run VBEM
for m = 2:200
    Estinfo = VBEMSmoother_full_nonlinear(spikeaccept,Constants,Estinfo,FieldMatrices,Basis,Thetainfo(m-1));
    Thetainfo(m) = VBM(Estinfo,Thetainfo(m-1),FieldMatrices,Constants,Basis,spikeaccept,Estmu,Estprec,Estkernel,Esttheta,Estrho);
%     Thetainfo(m).Meanprecmat = diag(Thetainfo(m).precisionest);
    save('LatestResults')
    %Breaking Condition
    if norm(Thetainfo(m).thetaest - Thetainfo(m-1).thetaest) < 0.005 ...
            && (abs(Thetainfo(m).muest - Thetainfo(m-1).muest) < 0.01) ...
            && max(Thetainfo(m-1).precisionest./Thetainfo(m).precisionest) < 1.05 ...
            && min(Thetainfo(m-1).precisionest./Thetainfo(m).precisionest) > 0.95 ...
            && norm(Thetainfo(m).kernelpar - Thetainfo(m-1).kernelpar) < 0.005
        break
    end
end

save(filename)
toc
%----------------------------------------------------------

function y = MultiIntForward2(x,Constants,Basis,Thetainfo)
beta = Constants.beta;
mu = Thetainfo.muest;
varmu = Thetainfo.varmu;
Uest = Basis.Basisvec*x;
y  = Constants.dt*Constants.ds^2*exp(mu+varmu/2)*beta^2*Basis.Basisvec'*(Basis.Basisvec.*repmat(exp(beta*Uest),1,Basis.nx));


function y = MultiIntForward(x,Constants,Basis,Thetainfo)

beta = Constants.beta;
mu = Thetainfo.muest;
varmu = Thetainfo.varmu;
Uest = Basis.Basisvec*x;
y =  Constants.dt*Constants.ds^2*exp(mu + varmu/2)*beta*Basis.Basisvec'*exp(beta*Uest);


