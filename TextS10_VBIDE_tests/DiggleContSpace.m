function [lambda] = DiggleContSpace(spikes,Constants)
% -------------------------------------------------------------
% Function Diggle
% inputs:
% outputs:
% Description: The Diggle estimator
%---------------------------------------------------------------

lambda = zeros(Constants.J,Constants.J);
frame = lambda;
buff = 1;
rabs = buff*Constants.ds;
% rabs = 3.2;
% buff = rabs/Constants.ds;
r = 4;
for i = 1:size(spikes,1)
    %Find x pixel
    [temp,x] = find((spikes(i,1) - Constants.s).^2 == min((spikes(i,1) - Constants.s).^2));
    %Find y pixel
    [temp,y] = find((spikes(i,2) - Constants.s).^2 == min((spikes(i,2) - Constants.s).^2));
    frame(y,x) = frame(y,x)+1;
end


for i = buff:Constants.J-buff
    for j = buff:Constants.J-buff
        lambda(j,i) = pixincircle(frame,i,j,r)/(pi*rabs^2)/Constants.dt;
    end
end
