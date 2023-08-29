%%% PAC-GS with compensation phase detected by AO
function [phase_para_noab,im]=WGS_para_corab_mod(System,mask,HStacks,phase_cor,Iw)
LZ=length(System.tilt_posi_x);
hologram = zeros(System.Nx,System.Ny);
for i = 1:LZ
    target = max(mask(:,:,i),0);
    hologram = hologram + ifft2(ifftshift(sqrt(target.*Iw(i)) .* exp(1i*randn(System.Nx,System.Ny))))./HStacks(:,:,i).*exp(1i*(phase_cor(:,:,i)));
end
im = System.source.*exp(1i * angle(hologram));
W = ones(LZ,1); 
for n = 1:System.maxiter
    
    index = 1:LZ;
    tempim = 0 * im;
    target_tem = zeros(LZ,1);
    phase_tem = zeros(System.Nx,System.Ny,LZ);

    for i = index
        imagez = fftshift(fft2(im.*exp(1i*(-phase_cor(:,:,i))).*HStacks(:,:,i)));
        I_mask = abs(imagez.*mask(:,:,i)).^2;
        target_tem(i) = sum(I_mask(:));
        phase_tem(:,:,i) = exp(1i * angle(imagez));
%         figure,imshow(I_mask,[])
    end
    target_tem = sum(target_tem(:))./(target_tem+System.GSoffset)/LZ;
    %     buf = (sum(total_spot_int(:)) / spot_number * intensity_weight(slice_ind) ./ (total_spot_int(:, slice_ind)) + 1e-3);
    W = W.*target_tem;
    for i = index
%         target = mask(:,:,i)*W(i) + System.GSoffset; % original WGS
         target = mask(:,:,i)*Iw(i) + System.GSoffset; % WGS by real intensity 
        imagez = sqrt(target) .* phase_tem(:,:,i);
        tempim = tempim + ifft2(ifftshift(imagez))./HStacks(:,:,i).*exp(1i*(phase_cor(:,:,i)));
    end
    im = System.source.*exp(1i * angle(tempim));
end
phase_para_noab = angle(im);
end