function [precision_x, precision_y, precision_z] = CRLB_cal(lambda,NA,Mag,pxl_size_cam,num_pxl,num_photon,bg_photon,n,z_range,z_interv,phase_mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                                     %
% This function calculates the theoretical loclaization precision of a certain point spread function with a set of input parameters.  %      
% It follows Poisson noise model discribed in Chao et al., J. Opt. Soc. Am. A 33, B36-B57 (2016).                                     %
%                                                                                                                                     %
% Input:	lambda			: wavelength in microns                                                                                   %
%			NA              : numerical aperture of objective                                                                         %
%			Mag             : magnification of imaging system                                                                         %
%			pxl_size_cam    : pixel size on the camera                                                                                %
%			num_pxl         : number of pixel in the field of view                                                                    %
%			num_photon      : expected number of sigmal photons per emitter                                                           %
%			bg_photon       : expected number of background photons per pixel                                                         %
%			n               : refractive index of the media                                                                           %
%			z_range         : range of interest along the optical axis                                                                %
%			z_interv        : Sampling interval along the optical axis                                                                %
%			phase_mask      : phase_mask, num_pxl x num_pxl complex numbered matrix				                                      %
%                             Note: effective area of the phase mask should match the pupil size                                      %
% Output:	zs				: Sampling points along the optical axis                                                                  %
%			precision_x		: Cramer-Rao lower bound in x (transverse direction)                                                      %
%			precision_y		: Cramer-Rao lower bound in y (transverse direction)                                                      %
%			precision_z		: Cramer-Rao lower bound in z (optical axis direction)                                                    %
% 																				                                                      %
%                                                                              														  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create the lateral coordinates
pxl_size 		= pxl_size_cam/Mag;														% effective pixel size, in microns
xx		 		= linspace(-pxl_size * num_pxl / 2, pxl_size * num_pxl / 2, num_pxl);	% lateral coordinates in the image plane, in microns

%% Create the Fourier plane coordinates
dx		 		= xx(2) - xx(1);  														% sampling period, in microns
fS		 		= 1 / dx;       														% spatial sampling frequency, inverse microns
df		 		= fS / num_pxl; 														% spacing between discrete frequency coordinates, inverse microns

[fx, fy]        = meshgrid(linspace(-df*(num_pxl/2-0.5), df*(num_pxl/2-0.5), num_pxl),...
				  linspace(-df*(num_pxl/2-0.5), df*(num_pxl/2-0.5), num_pxl)); 			% Fourier plane coordinates
[ftheta, fp]    = cart2pol(fx,fy);														% Fourier plane polar coordinates

%% Define the pupil aperture
fNA             = NA / lambda; 															% radius of the pupil, inverse microns
pupilAperture   = fp <= fNA; 															%Circular aperture
pupilsize		= sum(pupilAperture(end/2,:));


%% Calculate the 3D point spread function and Fisher information
PSFs			= [];            														% 3D PSF, initialize
PSFs_xplus		= [];																	% slightly shifted PSF in the direction of x
PSFs_yplus		= [];																	% slightly shifted PSF in the direction of y
PSFs_zplus		= [];																	% slightly shifted PSF in the direction of z
xplus			= 0.00001;																% shift in microns
yplus			= 0.00001;																% shift in microns
zplus			= 0.00001;																% shift in microns

bar				= waitbar(0);
zs				= -z_range:z_interv:z_range;
for ii = 1:length(zs)																	% loop to calculate the PSF at each depth
z 				= zs(ii);
DefocusPhase	= exp(1i.*2.*pi.*z.*sqrt((n./lambda).^2-fx.^2-fy.^2));					% phase term introduce by axial displacement
pupilFunc		= pupilAperture.*DefocusPhase.*phase_mask;								% pupil function
psf_a			= fftshift(fft2(pupilFunc));											% amplitude PSF
image			= abs(psf_a).^2;														% intensity PSF
image			= image./sum(image(:)).*num_photon; 									% PSF normalized and mutiplyed by number of photons
PSFs			= cat(3,PSFs,image); 													% stacking the PSF

% Calculate a shifted PSF along x
pupilFunc_xplus	= pupilAperture.*DefocusPhase.*phase_mask.*exp(1i.*2.*pi.*fx.*xplus);	% shift the psf by xplus microns
psf_a_xplus		= fftshift(fft2(pupilFunc_xplus));										% amplitude PSF
image_xplus		= abs(psf_a_xplus).^2;image_xplus=image_xplus./sum(image_xplus(:)).*num_photon; % PSF normalized and mutiplyed by number of photons
PSFs_xplus		= cat(3,PSFs_xplus,image_xplus);										% PSFs at different z 

% calculate a shifted PSF along y
pupilFunc_yplus	= pupilAperture.*DefocusPhase.*phase_mask.*exp(1i.*2.*pi.*fy.*yplus);	% shift the psf by xplus microns
psf_a_yplus		= fftshift(fft2(pupilFunc_yplus));										% amplitude PSF
image_yplus		= abs(psf_a_yplus).^2;image_yplus=image_yplus./sum(image_yplus(:)).*num_photon; % PSF normalized and mutiplyed by number of photons
PSFs_yplus		= cat(3,PSFs_yplus,image_yplus);										% PSFs at different z 

% calculate a shifted PSF along z
DefocusPhase_zplus	= exp(1i.*2.*pi.*(z+zplus).*sqrt((n./lambda).^2-fx.^2-fy.^2)); 
pupilFunc_zplus	= pupilAperture.*DefocusPhase_zplus.*phase_mask;						% pupil function
psf_a_zplus		= fftshift(fft2(pupilFunc_zplus));										% amplitude PSF
image_zplus		= abs(psf_a_zplus).^2;image_zplus=image_zplus./sum(image_zplus(:)).*num_photon;% PSF normalized and mutiplyed by number of photons
PSFs_zplus		= cat(3,PSFs_zplus,image_zplus);										% PSFs at different z 
waitbar(ii/length(zs),bar,'calculating sampling planes: ');
end
close(bar);
%% Partial derivative matrices for x, y, z
PPSFPx			= (PSFs_xplus-PSFs)./xplus;												% partial I_PSF partial_x
PPSFPy			= (PSFs_yplus-PSFs)./yplus;												% partial I_PSF partial_y
PPSFPz			= (PSFs_zplus-PSFs)./zplus;												% partial I_PSF partial_z

%% Calculate the Fisher information
Fisher_info_xx	= sum(sum(PPSFPx.^2./(PSFs+bg_photon),2),1);							% Fisher information diagnal element xx
Fisher_info_yy	= sum(sum(PPSFPy.^2./(PSFs+bg_photon),2),1);							% Fisher information diagnal element yy
Fisher_info_zz	= sum(sum(PPSFPz.^2./(PSFs+bg_photon),2),1);							% Fisher information diagnal element zz

%% Theoretical localization precision predicted by the Cramer-Rao lower bound
precision_x		= squeeze(squeeze(1./sqrt(Fisher_info_xx)));							% localization precision in x
precision_y		= squeeze(squeeze(1./sqrt(Fisher_info_yy)));							% localization precision in y
precision_z		= squeeze(squeeze(1./sqrt(Fisher_info_zz)));							% localization precision in z
																						

%% plot the precision
% figure;
% p				= plot((-size(precision_x,3)/2+1:size(precision_x,3)/2).*z_interv,squeeze(precision_x).*1000,...
% 				  (-size(precision_x,3)/2+1:size(precision_x,3)/2).*z_interv,squeeze(precision_y).*1000,'-.',...
% 				  (-size(precision_x,3)/2+1:size(precision_x,3)/2).*z_interv,squeeze(precision_z).*1000,'-');
% ylim([0,2.*1000.*max([max(squeeze(precision_x)),max(squeeze(precision_y)),max(squeeze(precision_z))])])
% p(1).LineWidth=2;
% p(1).Color='c';
% p(2).Color='m';
% p(2).LineWidth=2;
% p(3).LineWidth=2;
% p(3).Color='k';
% xlabel('z / um')
% ylabel('precision / nm')
% set(gca,'fontsize',14)
% xticks(-z_range:1:z_range)
% xlim([-z_range,z_range])

end