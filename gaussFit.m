function x = gaussFit(Target,Pos )
% This function fits a 2D Gaussian profile to the point image, returning parameters including Amplitude, x0, sigmax, y0, sigmay and angle
%% ----Area of interest----
MdataSize		= 18;																					% preset size of interest
Col				= round(Pos(1));																		% inital guess of the column number
Row				= round(Pos(2));																		% inital guess of the row number
if Row-MdataSize/2<1||Row+MdataSize/2>size(Target,1)||Col-MdataSize/2<1||Col+MdataSize/2>size(Target,2)
    x			= [0 0 0 0 0 0];
    return
end
Z				= Target(Row-MdataSize/2:Row+MdataSize/2,Col-MdataSize/2:Col+MdataSize/2);				% crop
%% Initial guess of parameters
x0				= [max(Z(:)),Pos(1),1,Pos(2),1,0];														% inital guess parameters [Amplitude, x0, sigmax, y0, sigmay, angle(in rad)]
%% Set coordinates
[X,Y]			= meshgrid(linspace(Col-MdataSize/2,Col+MdataSize/2,MdataSize+1),linspace(Row-MdataSize/2,Row+MdataSize/2,MdataSize+1));
xdata			= zeros(size(X,1),size(Y,2),2);
xdata(:,:,1)	= X;
xdata(:,:,2)	= Y;

%% --- Fit---------------------
% define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
lb				= [0,Col-MdataSize/2,0,Row-MdataSize/2,0,-pi/4];
ub				= [max(Z(:)),Col+MdataSize/2,MdataSize/2,Row+MdataSize/2,MdataSize/2,pi/4];%(MdataSize/2)^2
opts			= optimset('Display','off');
[x,resnorm,residual,exitflag] = lsqcurvefit(@Gauss_func2D,x0,xdata,Z,lb,ub,opts);
end

