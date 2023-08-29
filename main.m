%%% AO test based on Zernike mode;
%%% aberration compensation base on PAC-GS / add single or mean
%%% compensation phase with original GS hologram
%%% JC,2023/8/29

% framework:
% 1. original hologram of multiple foci input
% 2. AO test for each focus
% 3. PAC-GS and other compensation methods for beaded-ring patterns and
% foci compensation

% clear 
% clear all
% close all

load('ori_3pattern.mat');
load('base_Zernike_phase_upto27_fullpupil_withouttilt.mat');
%% 1. load orinigal hologram
sca
PsychDefaultSetup(2);
screens = Screen('Screens');
Screen('Preference','SkipSyncTests', 1);
screenNumber = 2; %SLM number
black = BlackIndex(screenNumber);
Screen('Preference', 'ConserveVRAM', 64);
Screen('Preference', 'SkipSyncTests', 1); % some internal self-tests and calibrations were skipped
[window, windowRect] = PsychImaging('OpenWindow',screenNumber,black); 
topPriorityLevel = MaxPriority(window);
ifi = Screen('GetFlipInterval', window)*2;
waitframes = 2;
Priority(topPriorityLevel);
vbl = Screen('Flip', window);
phase_noab = imread([loc,'focus_ori.bmp']);
Screen('PutImage', window, phase_noab);
Screen('Flip', window);

%% 2. sCMOS setting
imaqtool %Open Image Acquisition Toolbox, Determine where each focus is located, then close Image Acquisition Toolbox
wPosition = [415 389 1231 1326]; % cut the target plane to reduce memory usage, the numbers are read out in imaqtool, which is determined by the location of patterns
ROIPosition = [1086 435 95 83;1472 1494 83 87;500 1358 83 95]; % FOV for AO test of each focus
% close Image Acquisition Toolbox

vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');
src = getselectedsource(vid);
src.PixelType = 'mono16';
src.ExposureTime = 0.01;
vid.ROIPosition = wPosition; % whole FOV
vid.TriggerRepeat = Inf;
start(vid);
%% 3. AO test
tic
pattern_num = length(Depths);
phase_cor_AO = zeros(Nx,Nx,pattern_num); 
Zernike_ind = [1,2,9,3:8,10:24];
coef_scaling = [2,0.5,0.2];
iter = length(coef_scaling);
r=5; % area for signal accumulation
I = zeros(length(Zernike_ind),10,iter);
coef_best = zeros(length(Zernike_ind),iter,pattern_num);
frame = getsnapshot(vid);
thre_f = 65535; % set threshold to block abnormal bright pixel of sCMOS (sometimes need)
[lm,ln]=size(frame);
[sm,sn]=size(ROIPosition);
for mm = 1:sm
    phase_noab = imread([loc,'focus_ori.bmp']);
    Screen('PutImage', window, phase_noab);
    Screen('Flip', window);
    frame = zeros(lm,ln);
    for pp = 1:3
        frame_tem = getsnapshot(vid);
%         frame_tem(find(frame_tem>thre_f))=0; block abnormal bright pixel        
        frame = frame+double(frame_tem);
%          pause(0.05)
    end
    frame = frame/3;
    frame_roi = frame(ROIPosition(mm,2)-wPosition(2):ROIPosition(mm,2)-wPosition(2)+ROIPosition(mm,4),ROIPosition(mm,1)-wPosition(1):ROIPosition(mm,1)-wPosition(1)+ROIPosition(mm,3));
%     fprintf(['pattern',num2str(mm),',iteration',num2str(i),',Zernikr_ind',num2str(j),',value',num2str(max(frame_roi(:))),'\n']);
%     thre_f = 0.5*max(frame_roi(:));
    [BW]=target_mask(frame_roi,r);
%                 figure,imshow(frame_roi,[])   
    frame_eff = immultiply(frame_roi,BW);
    sig2 = frame_eff;
    Im_last = sum(sig2(:));
for i = 1:iter   
    for j = 1:length(Zernike_ind)
        coef = linspace(-coef_scaling(i),coef_scaling(i),10);
        coef_l = linspace(-coef_scaling(i),coef_scaling(i),100);
        I_test = zeros(1,length(coef));
        for k = 1:length(coef)     
            phase_tem = zeros(Nx,Ny);
            phase_Z = 2*pi*coef(k)*phase_Zer(:,:,Zernike_ind(j)); 
            phase_Z2 = phase_Z+phase_cor_AO(:,:,mm)+phase_focus;
            phase_Z_wrap = mod(phase_Z2,2*pi);
            phase_tem(:,Ny/2-Nx/2:Ny/2+Nx/2-1) = phase_Z_wrap;
            phase_00 = double(uint8((phase_tem)/ 2 / pi * 255))+corrected_pat;
            phase_01 = mod(phase_00, 256);
            phase_02 = uint8(phase_01);
            Screen('PutImage', window, phase_02);
            Screen('Flip', window);
            pause(0.08)            
            frame = zeros(lm,ln);
            for pp = 1:3
                frame_tem = getsnapshot(vid);
%                 frame_tem(find(frame_tem>thre_f))=0; 
                frame = frame+double(frame_tem);
%                 pause(0.05)
            end
            frame = frame/3;           
            frame_roi = frame(ROIPosition(mm,2)-wPosition(2):ROIPosition(mm,2)-wPosition(2)+ROIPosition(mm,4),ROIPosition(mm,1)-wPosition(1):ROIPosition(mm,1)-wPosition(1)+ROIPosition(mm,3));
%             fprintf(['pattern',num2str(mm),',iteration',num2str(i),',Zernikr_ind',num2str(j),',value',num2str(max(frame_roi(:))),'\n']);
            [BW]=target_mask(frame_roi,r);
%             figure,imshow(frame_roi,[])
            frame_eff = immultiply(frame_roi,BW);
            sig2 = frame_eff;
            I(j,k,i) = sum(sig2(:));
            I_test(k) = sum(sig2(:));
        end
        [Im,ind] = max(I_test);
        if Im>Im_last
        coef_best_tem = coef(ind);
        coef_best(j,i,mm) = coef_best_tem;
        phase_cor_AO(:,:,mm) = phase_cor_AO(:,:,mm)+2*pi*coef_best_tem*phase_Zer(:,:,Zernike_ind(j));
        Im_last = Im;
        end        
    end
end
end

%% 4. Phase compensation by PAC-GS and phase adding methods
% PAC-GS
[phase_para_corab_f_AO]=WGS_para_corab_mod(System,I_target_focus,HStacks,phase_cor_AO,Iw_focus);
[phase_para_corab_BR_AO]=WGS_para_corab_mod(System,I_target_BR,HStacks,phase_cor_AO,Iw_BR);
toc
% image acquisition
% after AO
phase_Z_wrap = mod(phase_para_corab_f_AO,2*pi);
phase_tem(:,Ny/2-Nx/2:Ny/2+Nx/2-1) = phase_Z_wrap;
phase_00 = double(uint8((phase_tem)/ 2 / pi * 255)) + corrected_pat;
phase_01 = mod(phase_00, 256)*wave_ratio;
phase_02 = uint8(phase_01);
Screen('PutImage', window, phase_02);
Screen('Flip', window);
pause(0.5)
frame_a = getsnapshot(vid);
PAC_GS_focus = max(frame_a(:))
% figure,imshow(frame_a,[])
saveastiff(frame_a,'after_AO_focus.tif');
phase_Z_wrap = mod(phase_para_corab_BR_AO,2*pi);
phase_tem(:,Ny/2-Nx/2:Ny/2+Nx/2-1) = phase_Z_wrap;
phase_00 = double(uint8((phase_tem)/ 2 / pi * 255)) + corrected_pat;
phase_01 = mod(phase_00, 256)*wave_ratio;
phase_02 = uint8(phase_01);
Screen('PutImage', window, phase_02);
Screen('Flip', window);
pause(0.5)
frame_a = getsnapshot(vid);
PAC_GS_BR = max(frame_a(:))
saveastiff(frame_a,'after_AO_BR.tif');

% before AO
phase_Z_wrap = mod(phase_focus,2*pi);
phase_tem(:,Ny/2-Nx/2:Ny/2+Nx/2-1) = phase_Z_wrap;
phase_00 = double(uint8((phase_tem)/ 2 / pi * 255)) + corrected_pat;
phase_01 = mod(phase_00, 256)*wave_ratio;
phase_02 = uint8(phase_01);
Screen('PutImage', window, phase_02);
Screen('Flip', window);
pause(0.5)
frame_b = getsnapshot(vid);
ab_focus = max(frame_b(:))
% figure,imshow(frame_b,[])
saveastiff(frame_b,'before_AO_focus.tif');
phase_Z_wrap = mod(phase_BR,2*pi);
phase_tem(:,Ny/2-Nx/2:Ny/2+Nx/2-1) = phase_Z_wrap;
phase_00 = double(uint8((phase_tem)/ 2 / pi * 255)) + corrected_pat;
phase_01 = mod(phase_00, 256)*wave_ratio;
phase_02 = uint8(phase_01);
Screen('PutImage', window, phase_02);
Screen('Flip', window);
pause(0.5)
frame_b = getsnapshot(vid);
ab_BR = max(frame_b(:))
% figure,imshow(frame_b,[])
saveastiff(frame_b,'before_AO_BR.tif');

% GS+mean compensation phase
phase_cor_AO_mean = sum(phase_cor_AO,pattern_num)/pattern_num;
phase_Z_wrap = mod(phase_focus+phase_cor_AO_mean,2*pi);
phase_tem(:,Ny/2-Nx/2:Ny/2+Nx/2-1) = phase_Z_wrap;
phase_00 = double(uint8((phase_tem)/ 2 / pi * 255)) + corrected_pat;
phase_01 = mod(phase_00, 256)*wave_ratio;
phase_02 = uint8(phase_01);
Screen('PutImage', window, phase_02);
Screen('Flip', window);
pause(0.5)
frame_m = getsnapshot(vid);
meancor_focus = max(frame_m(:))
% figure,imshow(frame_m,[])
saveastiff(frame_m,'meancor_AO_focus.tif');
phase_cor_AO_mean = sum(phase_cor_AO,pattern_num)/pattern_num;
phase_Z_wrap = mod(phase_BR+phase_cor_AO_mean,2*pi);
phase_tem(:,Ny/2-Nx/2:Ny/2+Nx/2-1) = phase_Z_wrap;
phase_00 = double(uint8((phase_tem)/ 2 / pi * 255)) + corrected_pat;
phase_01 = mod(phase_00, 256)*wave_ratio;
phase_02 = uint8(phase_01);
Screen('PutImage', window, phase_02);
Screen('Flip', window);
pause(0.5)
frame_m = getsnapshot(vid);
meancor_BR = max(frame_m(:))
% figure,imshow(frame_m,[])
saveastiff(frame_m,'meancor_AO_BR.tif');

% GS+single compensation phase
for kk = 1:pattern_num
phase_Z_wrap = mod(phase_focus+phase_cor_AO(:,:,kk),2*pi);
phase_tem(:,Ny/2-Nx/2:Ny/2+Nx/2-1) = phase_Z_wrap;
phase_00 = double(uint8((phase_tem)/ 2 / pi * 255)) + corrected_pat;
phase_01 = mod(phase_00, 256)*wave_ratio;
phase_02 = uint8(phase_01);
Screen('PutImage', window, phase_02);
Screen('Flip', window);
pause(0.5)
frame_b = getsnapshot(vid);
max(frame_b(:))
% figure,imshow(frame_b,[])
saveastiff(frame_b,['cor_AO_focus',num2str(kk),'.tif']);

phase_Z_wrap = mod(phase_BR+phase_cor_AO(:,:,kk),2*pi);
phase_tem(:,Ny/2-Nx/2:Ny/2+Nx/2-1) = phase_Z_wrap;
phase_00 = double(uint8((phase_tem)/ 2 / pi * 255)) + corrected_pat;
phase_01 = mod(phase_00, 256)*wave_ratio;
phase_02 = uint8(phase_01);
Screen('PutImage', window, phase_02);
Screen('Flip', window);
pause(0.5)
frame_b = getsnapshot(vid);
max(frame_b(:))
% figure,imshow(frame_b,[])
saveastiff(frame_b,['cor_AO_BR',num2str(kk),'.tif']);
end
save([loc,'AO_3 pattern.mat']);
delete(vid);