function y = gaussinner(s1,s2,phi)
% s1, s2 are generated from meshgrid of the space discretization
% phi, the Localized kernel for each kernel

    nx = size(phi,3);   % the number of kernels
    y = zeros(nx,nx);   % init the relation between kernels?
    
    for i = 1:nx
        for j = 1:nx
            % trapz is integration function
            %
            % TODO, quite difficult to understand here
            % phi(:,:,i).*phi(:,:,j), len(s1)*len(s2) matrix
            %
            %
            y(i,j) = trapz( s2(:,1),(trapz(s1(1,:),phi(:,:,i).*phi(:,:,j),2)) );
        end
    end
    
end