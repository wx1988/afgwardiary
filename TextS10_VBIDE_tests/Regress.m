
function y = Regress(lambda,Basis,mu,beta)

phivec = (reshape(Basis.phi,[],Basis.nx))';
y = zeros(Basis.nx,size(lambda,3));
lambda = lambda + 0.001;
for i = 1:size(lambda,3)
    lambdavec = reshape(lambda(:,:,i),[],1);
    y(:,i) = inv(phivec*phivec')*phivec*((log(lambdavec) - mu)/beta);
end