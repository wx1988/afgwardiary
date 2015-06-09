function y = SingleIntForward(x,Constants,Basis,Thetainfo)
beta = Constants.beta;
mu = Thetainfo.muest;
varmu = Thetainfo.varmu;
Uest = Basis.Basisvec*x;
y = Constants.dt*Constants.ds^2*exp(mu + varmu/2)*sum(exp(beta*Uest));
