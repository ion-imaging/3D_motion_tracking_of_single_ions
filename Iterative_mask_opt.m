%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                                     %
% This script performs an iterative optimization for the helical phase mask. It starts from the Fresnel-zone-based method described   % 
% in Berlich et al., Opt. Express 26, 4873-4891 (2018) and is optimized using an adapted version of the Gerchberg-Saxton algorithm.   %
% Note that the initial phase mask defined on the Fresnel zones can be directly employed in experiment on an spatial light modulator, %
% however, this optimization process reduces the the number of diffraction edges in the mask, making it more suitable for fabricated  %
% masks on glass substrates; besides, the optimization process can also boost the two PSF lobes at different rotation angles.         %
%                                                                                                                                     %
%                                                                              														  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;clear all;
%% Setup the simulation parameters
lambda				= 0.397;																% wavelength in microns
NA        			= 0.4;																	% numerical aperture of the objective
Mag         		= 40;																	% magnification of the imaging system
pxl_size 			= 6.5./Mag;																% effective pixel size inmicrons
num_pxl 			= 2000;																	% number of pixels in the camera
num_photon 			= 1000;																	% number of signal photons 
bg_photon  			= 5;																	% background photon counts per pixel  
n         			= 1;																	% refractive index of the media

z_range				= 4.8;																	% depth range of interest in microns
z_interv			= 1.2;																	% interval between sampling planes in microns
zs					= -z_range:z_interv:z_range;											% sampling planes along the optical axis

num_iter			= 10;																	% number of iterations

%% Create the image plane coordinates
xx					= linspace(-pxl_size*(num_pxl/2-0.5), pxl_size*(num_pxl/2-0.5), num_pxl);% in microns
x_pxl               = linspace(1,num_pxl,num_pxl);                                          % 
[x_pxl,y_pxl]       = meshgrid(x_pxl,x_pxl);                                                % xy coordinates in pixels

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
phase_mask			= exp(1i.*DH_phase_Fresnl(fp,ftheta,fNA,4,2,0.5));						% example initial phase mask
figure; imshow(angle(phase_mask),[]);														% show the phase mask

%% Iterative optimization
for kk = 1:num_iter
pupilFunc_avr		= zeros(num_pxl,num_pxl);												% averaged pupil function

for zi = 1:length(zs)
    z				= zs(zi);
    DefocusPhase	= exp(1i.*2.*pi.*z.*sqrt((n./lambda).^2-fx.^2-fy.^2)); 					% phase term introduced by axial displacement of the emitter
    pupilFunc		= pupilAperture.*DefocusPhase.*phase_mask.*(fx>0);						% pupil function left half
    psf_a			= fftshift(fft2(pupilFunc)).*dx;										% amplitude PSF
    image			= abs(psf_a).^2;														% intensity PSF
	image			= image./sum(image(:)).*num_photon;										% PSF normalized and mutiplyed by number of photons
    %imshow(image,[]);
    
    [Mag,I]			= max(image(:));
    [maxY, maxX]	= ind2sub(size(image),I);												% locate the peak of the PSF lobe
    fit_res			= gaussFit(image,[maxX, maxY]);											% fit to 2D Gaussian function
    gauss_filter	= Gauss_func2D([fit_res(1),fit_res(2),fit_res(3),fit_res(4),fit_res(5),fit_res(6)],cat(3,x_pxl,y_pxl)); % use the fitting result to construct a 2D filter
    gauss_filter	= gauss_filter./max(gauss_filter(:));									% normalize the filter
    psf_a2			= abs(psf_a).*gauss_filter.^0.5.*exp(1i.*angle(psf_a));					% apply the filter to the amplitude PSF
    
    %psf_a2 = abs(psf_a).*exp(1i.*angle(psf_a));
    pupilFunc		= ifft2(ifftshift(psf_a2))./dx./DefocusPhase;							% inverse Fourier transform to get the pupil function
    
    %figure;imshow(pupilAperture.*angle(pupilFunc),[]);%.*(fx<0)

    pupilFunc_avr	= pupilFunc_avr + pupilFunc;											% sum for all corresponding sampling planes
end

pupilFunc_avr		= pupilFunc_avr./length(zs);											% averaged
	
phase_mask			= exp(1i.*angle(pupilFunc_avr));										% phase mask for next iteration

imshow(angle(phase_mask).*pupilAperture.*(fx>0)+imrotate(angle(phase_mask).*pupilAperture.*(fx>=0),180),[]); % construct the two halves and show the current phase mask
%imshow(image,[]);

title(['Iteration ', num2str(kk), ' in ', num2str(num_iter), '.'],'FontSize',18);pause(0.1);
end

phase_0				= angle(phase_mask).*pupilAperture.*(fx>0)+imrotate(angle(phase_mask).*pupilAperture.*(fx>=0),180); % optimized phase mask (wrapped)

%% Phase unwrapping
phase_1				= phase_0;
for kk= 1:num_pxl
    phase_1(:,kk)	= unwrap(phase_0(:,kk));												% 1D phase unwrapping
    %kk
end
%imshow(phase_1,[])

phase_2				= phase_1;
for kk= 1:num_pxl
    phase_2(:,kk)	= phase_1(:,kk)-phase_1(uint16(round(num_pxl/2))+2,kk);										% unwrap the other dimension
    %kk
end
phase_3				= phase_2.*pupilAperture;
%imshow(phase_3,[])
%save('phase_3_0','phase_3_0');

%% Adding linear phase ramp to seperate unmodulated light (for phase masks with non-ideal fill factors)
phase_mask			= exp(1i.*(phase_3 + 3.5.*fy));											
imshow(angle(phase_mask).*pupilAperture,[]);													% show the final optimized phase mask
 