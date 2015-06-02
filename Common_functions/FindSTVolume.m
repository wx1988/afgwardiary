function Volume = FindSTVolume(s1,s2,t,GaussKernel,mu)
% FINDSTVOLUME, Find the volume (numerically) of a 3D GRBF centred on mu on the domain (s1,s2,t) 
% 
% See also MEAN, DIFF, MESHGRID, TRAPZ
%
% 

Volume = 0;
% The diff function will calculate the difference between two items. 
dt = mean(diff(t));

% TODO, the ds1 and ds2 are not used later
ds1 = mean(diff(s1));
ds2 = mean(diff(s2));

% meshgrid will generate two matrix of len(s2) * len(s1)
% the first matrix with each row the same as s1
% the second matrix with each column the same as s2
[S1,S2] = meshgrid(s1,s2);

% TODO, cannot understand here
for i = 1:length(t)
    % The trapz is trying to calculate the integral
    Volume = Volume + trapz(s2,trapz(s1,GaussKernel(S1,S2,t(i),mu(1),mu(2),mu(3))))*dt;
end