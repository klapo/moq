function [ym,ys,yn] = bindata2(y,x1,x2,x1rg,x2rg)
%function [ym,ys,ns] = bindata2(y,x1,x2,x1rg,x2rg)
%Computes:
%ym(ii,jj) = mean(y(x1>=x1rg(ii) & x1 < x1rg(ii+1) & x2>=x2rg(jj) & x2 < x2rg(jj+1))
%for every ii, jj
%If a bin is empty it returns nan for that bin
%using a fast algorithm which uses no looping
%Also returns ys, the standard deviation of y using binning (useful for r^2
%calculations). Example:
%
% x = randn(500,2);
% y = sum(x.^2,2) + randn(500,1);
% xrg = linspace(-3,3,10)';
% [ym,yb] = bindata2(y,x(:,1),x(:,2),xrg,xrg);
% subplot(1,2,1);plot3(x(:,1),x(:,2),y,'.');
% subplot(1,2,2);h = imagesc(xrg,xrg,ym);
% set(h,'AlphaData',~isnan(ym)); box off;
%
%Based on function written by Patrick Mineault
%Refs: http://xcorr.net/?p=3326
%      http://www-pord.ucsd.edu/~matlab/bin.htm
[~,whichedge1] = histc(x1,x1rg(:)');
[~,whichedge2] = histc(x2,x2rg(:)');

bins1 = min(max(whichedge1,1),length(x1rg)-1);
bins2 = min(max(whichedge2,1),length(x2rg)-1);

bins = (bins2-1)*(length(x1rg)-1)+bins1;


%% Mean
xpos = ones(size(bins,1),1);
ns = sparse(bins,xpos,1,(length(x1rg)-1)*(length(x2rg)-1),1);
ysum = sparse(bins,xpos,y,(length(x1rg)-1)*(length(x2rg)-1),1);
ym = full(ysum)./(full(ns));
%% Bin count
yn = full(ns);
%% Standard Deviation
ydiff = y-(ysum(bins(:),1)./ns(bins(:),1));
ydiff = sparse(bins,xpos,ydiff.^2,(length(x1rg)-1)*(length(x2rg)-1),1);
ys = sqrt(full(ydiff)./(full(ns)-1));
ys(ys==0) = NaN;

%% Output
ym = reshape(ym,length(x1rg)-1,length(x2rg)-1);
ys = reshape(ys,length(x1rg)-1,length(x2rg)-1);
yn = reshape(yn,length(x1rg)-1,length(x2rg)-1);
end
