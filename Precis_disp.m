%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                                     %
% This script performs a comparison between double helical PSF and the conventional Airy disk method by calculating their theoretical %
% localization precision. Default parameters are set to meet our experimental conditions for the ion tracking study. 				  %
%                                                                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
%% Setup the simulation parameters
lambda				= 0.397;																% wavelength in microns
NA        			= 0.4;																	% numerical aperture of the objective
Mag         		= 36;																	% magnification of the imaging system
pxl_size_cam 		= 13;                                                          			% effective pixel size inmicrons
num_pxl 			= 1000;																	% number of pixels in the camera
num_photon 			= 22000;																% number of signal photons 
bg_photon  			= 10;																	% background photon counts per pixel  
n         			= 1;																	% refractive index of the media

%% sampling on the depth range of interest
vRange      = 10;      % z range in micron
vInterv     = 0.2;   % z interval in micron
v_sampl = -vRange:vInterv:vRange;

%% Create the image plane coordinates
pxl_size = pxl_size_cam/Mag;                                                                % effective pixel size
xx					= linspace(-pxl_size*(num_pxl/2-0.5), pxl_size*(num_pxl/2-0.5), num_pxl);% in microns
%% Create the Fourier plane coordinates
dx					= xx(2) - xx(1); 														% sampling period, microns
fS					= 1 / dx;        														% spatial sampling frequency, inverse microns
df					= fS / num_pxl;															% spacing between discrete frequency coordinates, inverse microns
[fx, fy]			= meshgrid(linspace(-df*(num_pxl/2-0.5), df*(num_pxl/2-0.5), num_pxl),...
					  linspace(-df*(num_pxl/2-0.5), df*(num_pxl/2-0.5), num_pxl)); 			% Fourier plane coordinates
[ftheta, fp]		= cart2pol(fx,fy);														% Fourier plane polar coordinates
%% Define the pupil aperture
fNA          		= NA / lambda;															% radius of the pupil, inverse microns
pupilAperture		= fp <= fNA;															% circular aperture
pupilsize			= sum(pupilAperture(end/2,:));
%% Define initial phase mask
DHphase_mask		= exp(1i.*DH_phase_Fresnl(fp,ftheta,fNA,4,2,0.9));						% example initial phase mask
%figure; imshow(angle(DHphase_mask),[]);													% show the phase mask

%% Calculate the theoretical precision in h(transverse), z(transverse), v(optical axis) directions
[DHprecision_z, DHprecision_h, DHprecision_v] = CRLB_cal(lambda,NA,Mag,pxl_size_cam,num_pxl,num_photon,bg_photon,n,vRange,vInterv,DHphase_mask); % precision for double helical PSF
[ADprecision_z, ADprecision_h, ADprecision_v] = CRLB_cal(lambda,NA,Mag,pxl_size_cam,num_pxl,num_photon,bg_photon,n,vRange,vInterv,pupilAperture); % precision for conventional Airy disk

%% Plot the precision curves
figure;
p=plot(v_sampl,DHprecision_z.*1000,'-',v_sampl,DHprecision_h.*1000,'-',v_sampl,DHprecision_v.*1000,'-');hold on
ylim([0,1.5.*1000.*max(DHprecision_v)]);
p(1).LineWidth=2;
p(1).Color='c';
p(2).Color='m';
p(2).LineWidth=2;
p(3).LineWidth=2;
p(3).Color='k';

p=plot(v_sampl,ADprecision_z.*1000,'-.',v_sampl,ADprecision_h.*1000,'--',v_sampl,ADprecision_v.*1000,'-');
ylim([0,1.5.*1000.*max(DHprecision_v)]);
p(1).LineWidth=2;
p(1).Color='g';
p(2).Color='r';
p(2).LineWidth=2;
p(3).LineWidth=2;
p(3).Color='b';

title('theoretical localization precision')
xlabel('z / um')
ylabel('precision / nm')
set(gca,'fontsize',14)
xlim([-vRange,vRange])
grid on
legend('helical z','helical h', 'helical v','conventional z','conventional h','conventional v');