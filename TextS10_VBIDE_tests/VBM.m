
function [Thetainfo] = VBM(Estinfo,Thetainfo,FieldMatrices,Constants,Basis,spikes,Estmu,Estprec,Estkernel,Esttheta,Estrho)

PSIx = FieldMatrices.PSIx;
PSIxinv = inv(PSIx);
SigmaW = FieldMatrices.W;
d = Constants.d;
dt = Constants.dt;
ds = Constants.ds;
N = Constants.N;
PSIxinvV = PSIxinv*FieldMatrices.V(:,:,1);
rho =  Thetainfo(1).rho;
varrho = Thetainfo(1).varrho;
Aest = rho*Thetainfo.kernelpar*PSIxinvV;

%Estimate precision
alpha = 5;
beta = 0.2;

precest = repmat(1/(0.2^2),Basis.nx,1);
varprec = repmat(0,Basis.nx,1);
if Estprec == 'y'
        for j = 1:Basis.nx
            k1 = 0;
            for i = 3:N-2
                k1 = k1+ Estinfo(i).W(j,j) + rho^2*Estinfo(i-1).W(j,j) + dt^2*(Thetainfo(1).vartheta(j,j) + Thetainfo(1).thetaest(j)^2) ...
                    - 2*rho*Estinfo(i-1).S(j,j) + 2*dt*Thetainfo(1).thetaest(j)*Estinfo(i-1).xestRTS(j)*rho ...
                    - 2*dt*Thetainfo(1).thetaest(j)*Estinfo(i).xestRTS(j);
            end
            alpha0 = (Constants.N-4)/2;%alpha0 = alpha + (Constants.N-4)/2;
            beta0 = k1/2;%beta0 = beta + k1/2;
            precest(j) = alpha0/beta0;
            varprec(j) = alpha0/beta0^2;
        end
        Meanprecmat = diag(precest);
    
    %Whole Wishart distribution
    dofprior = 10;
    Prevmatprior = 2/dofprior*eye(Basis.nx);
    Gamma = 0;
    for i = 3:N-2
        Gamma = Gamma + Estinfo(i).W + (rho^2 + varrho)*(Thetainfo(1).kernelpar^2 + Thetainfo(1).varkernelpar)*PSIxinvV*Estinfo(i-1).W*PSIxinvV' ...
            + Thetainfo(1).vartheta + Thetainfo(1).thetaest*Thetainfo(1).thetaest' ...
            - Aest*Estinfo(i-1).S - Estinfo(i-1).S'*Aest' ...
            - Estinfo(i).xestRTS*Thetainfo(1).thetaest' - Thetainfo(1).thetaest* Estinfo(i).xestRTS'  ...
            + Aest*Estinfo(i-1).xestRTS*Thetainfo(1).thetaest' + Thetainfo(1).thetaest*Estinfo(i-1).xestRTS'*Aest';
    end
    Precmatpost = inv(inv(Prevmatprior) + Gamma);
    Meanprecmat = (dofprior + N-4)*Precmatpost;
else
    Meanprecmat = Thetainfo(1).Meanprecmat;
end



%Estimate theta
if Esttheta == 'y'
    v = zeros(d,1);
    for i = 3:Constants.N-2
        v(:,1) = v(:,1) + dt*inv(SigmaW)*(Estinfo(i).xestRTS - Aest*Estinfo(i-1).xestRTS);
    end
    vartheta = SigmaW./dt^2/(Constants.N-4);
    thetaest = vartheta*v;
else
   thetaest = Thetainfo.thetaest;
   vartheta = Thetainfo.vartheta;
end



%Estimaterho
if Estrho == 'y'
    rhoprior = 1;
    varrhoprior = 1;
    v = 0;
    Upsilon = 0;
    for i = 3:N-2
        v = v + trace(Thetainfo(1).Meanprecmat*Estinfo(i-1).S' - Thetainfo(1).Meanprecmat*Thetainfo.thetaest*Estinfo(i-1).xestRTS');
        Upsilon = Upsilon + trace(Thetainfo(1).Meanprecmat*Estinfo(i-1).W);
    end
    Thetainfo(1).varrho = inv(inv(varrhoprior) + Upsilon);
    Thetainfo(1).rho = Thetainfo(1).varrho*(inv(varrhoprior)*rhoprior + v);
end


%EstimateKernel
if Estkernel == 'y'
    kernelprior = 1;
    varkernelprior = 1;
    v = 0;
    Upsilon = 0;
    for i = 3:N-2
        v = v + trace(Thetainfo.Meanprecmat*PSIxinv*FieldMatrices.V(:,:,1)*Estinfo(i).S' - Thetainfo.Meanprecmat*Thetainfo.thetaest*Estinfo(i).xestRTS'*FieldMatrices.V(:,:,1)'*PSIxinv);
        Upsilon = Upsilon + trace(FieldMatrices.V(:,:,1)'*PSIxinv*Thetainfo.Meanprecmat*PSIxinv*FieldMatrices.V(:,:,1)*Estinfo(i).W);
    end
    Thetainfo(1).varkernelpar = inv(inv(varkernelprior) + Upsilon);
    Thetainfo(1).kernelpar = Thetainfo(1).varkernelpar*(inv(varkernelprior)*kernelprior + v);
end

%Estimate mu
if Estmu == 'y'
    muprior = 0;
    Sigmamuprior = 10;
    totalspikenum = 0;
    totalvoidsum = 0;
    phitemp = zeros(Constants.J,Constants.J,Basis.nx);
    Basis_vec_mat = reshape(Basis.phi,[],Basis.nx);
    for i = 3:N-2
        Umean = sum(multiprod(Basis.phi,reshape(Estinfo(i).xestRTS,1,1,Basis.nx)),3);
        MySigma = Estinfo(i).PRTS; %CHANGE TO FULL COVARIANCE
        phitemp_vec = Basis_vec_mat*MySigma;
        phitemp = reshape(phitemp_vec,Constants.J,Constants.J,Basis.nx);
        %         for j = 1:Basis.nx
        %             phitemp(:,:,j) = sum(multiprod(Basis.phi,reshape(MySigma(:,j),1,1,Basis.nx)),3);
        %         end
        phisigma = sum(phitemp.*Basis.phi,3);
        totalvoidsum = totalvoidsum + dt*Constants.ds^2*sum(sum(exp(Constants.beta*Umean + Constants.beta^2*phisigma./2)));
        totalspikenum = totalspikenum + size(spikes(i).Coords,1);
    end
    fmu = @(mu) -(-(mu - muprior)^2/(2*Sigmamuprior) + mu*totalspikenum - exp(mu)*totalvoidsum);
    gradfmu = @(mu) -(-mu/Sigmamuprior + totalspikenum - exp(mu)*totalvoidsum);
    options = foptions;
    options(9) = 1;
    Thetainfo(1).muest = scg(fmu,muprior,options,gradfmu)';
    Thetainfo(1).varmu = 1/(1/Sigmamuprior + exp(Thetainfo(1).muest)*totalvoidsum);
end

Thetainfo(1).precisionest = precest;
Thetainfo(1).varprecision = varprec;
Thetainfo(1).Meanprecmat = Meanprecmat;
Thetainfo.thetaest = thetaest;
Thetainfo.vartheta = vartheta;