function gauss_2D = Gauss_func2D(x,xdata)
% This function takes the x coordinate xdata(:,:,1), y coordinate xdata(:,:,2) and parameters of x = [Amp, x0, wx, y0, wy, fi]; it returns a 2D Gaussian profile 
% 
xcoord			= xdata(:,:,1);
ycoord			= xdata(:,:,2);
xcoord_rot		= xcoord.*cos(x(6)) - ycoord.*sin(x(6));
ycoord_rot		= xcoord.*sin(x(6)) + ycoord.*cos(x(6));
x0rot 			= x(2).*cos(x(6)) - x(4).*sin(x(6));
y0rot 			= x(2).*sin(x(6)) + x(4).*cos(x(6));

gauss_2D		= x(1).*exp(   -((xcoord_rot-x0rot).^2./(2.*x(3)^2) + (ycoord_rot-y0rot).^2/(2.*x(5).^2) ));