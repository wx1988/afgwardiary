function y = pixincircle(A,xc,yc,radius)

% Engine
y=1:size(A,1);
x=1:size(A,2);
[X Y]=meshgrid(x,y);

% This assume the circle falls *entirely* inside the image
R2 = (X-xc).^2+(Y-yc).^2;
% c = contourc(x,y,R2,[0 0]+radius^2);
[c1,c2] = find(R2 < radius^2);
% c = round(c(:,2:end)); % pixels located ~ on circle
c(1,:) = c1;
c(2,:) = c2;
c = round(c(:,2:end)); % pixels located in circle
Ac = A(sub2ind(size(A),c(1,:),c(2,:))); % extract value
y = sum(Ac);