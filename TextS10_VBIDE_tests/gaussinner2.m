function y = gaussinner2(s1,s2,phi1,phi2)
nx = size(phi1,3);
y = zeros(nx,nx);
for i = 1:nx
    for j = 1:nx
        y(i,j) = trapz(s2(:,1),(trapz(s1(1,:),phi1(:,:,i).*phi2(:,:,j),2)));
    end
end