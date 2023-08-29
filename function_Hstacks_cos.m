function [ HStacks] = function_Hstacks_cos( System,z )
%Precompute Fresnel propagation kernels at depths z 

psXHolograph = System.lambda * System.focal_SLM * System.focal_obj / System.psSLM  / System.Nx / System.focal_tube;      % Pixel Size (resolution) at the scattered 3D region
psYHolograph = System.lambda * System.focal_SLM * System.focal_obj / System.psSLM  / System.Ny / System.focal_tube;      % Pixel Size (resolution) at the scattered 3D region

if System.useGPU ==1
    HStacks = zeros(System.Nx, System.Ny, numel(z), 'gpuArray');
else
    HStacks = zeros(System.Nx, System.Ny, numel(z));                     
end
for i = 1 : numel(z)
    [HStacks(:,:,i)] = function_GenerateFresnelPropagationStack_cos(System.Nx, System.Ny, z(i), System.lambda, System.useGPU, System.ObjNA, System.ObjRI);
end
end

