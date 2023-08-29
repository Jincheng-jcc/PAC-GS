function [HStack]=function_GenerateFresnelPropagationStack_cos(Nx,Ny,z, lambda, useGPU,NA,RI)
% Lambda, ps, z has unit meter.
% z is the distance from focal plane.
% lambda is wavelength
% ps is pixel size.
cx=(1:Nx) - (floor(Nx/2)+1);
cy=(1:Ny) - (floor(Ny/2)+1);

if useGPU
    cx = gpuArray(cx); cy = gpuArray(cy);
end
k0 = 2*pi/lambda;
[kx, ky] = meshgrid(cx, cy);
circle = (kx).^2+(ky).^2; 
R=Ny/2;
k = sqrt(kx.^2+ky.^2);
dkxy = NA/RI/floor(Nx/2);
sinthetaObj = k.*dkxy;
sinthetaObj(find(circle>R^2)) = 1;
costhetaObj = eps+sqrt(1-(sinthetaObj.^2));
phid = (sqrt(-1)*k0*RI).*costhetaObj;
HStack = exp(phid * z);

end