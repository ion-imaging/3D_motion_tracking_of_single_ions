function DHphase = DH_phase_Fresnl(rho0,theta,Rpup,L,N,ep)
% this function generates the double helical phase mask defined on circular
% Fresnel zones [Berlich et al., Opt. Express 26, 4873-4891 (2018)]
rho = rho0./Rpup;
DHphase = zeros(size(rho));
for i=1:L
    DHphase(rho>=(i-1)^ep/L^ep & rho<i^ep/L^ep) = ((i-1).*N+1).*theta(rho>=(i-1)^ep/L^ep & rho<i^ep/L^ep);
end
% A=DHphase;
% DHphase = angle(exp(1i.*DHphase));
% figure;imshow(A,[]);
% figure;imshow(DHphase,[]);
end

