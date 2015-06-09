
function V = FindV(s1,s2,Basis,KernelBasis)
nx = Basis.nx;
d = KernelBasis.d;
V = zeros(d,nx,nx);
for i = 1:d
    for j = 1:nx
        for k = 1:nx
            V(j,k,i) = trapz(s2(:,1),(trapz(s1(1,:),Basis.phi(:,:,j).*KernelBasis.PHI(:,:,i,k),2)));
        end
    end
end