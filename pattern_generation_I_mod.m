%%% wGS for original multiple foci and beaded-ring pattern generation
%%% JC, 23/8/29
%% System parameter
close all
clear all
System.verbose=1;           % 1 or 0    Set this value to 1 to display activity, 0 otherwise
System.lambda = 1.04e-6;    % unit:m
System.psSLM = 20e-6;       % meters    SLM pixel dimensions
System.Nx = 600;            % int       Number of pixels in X direction
System.Ny = 600;            % int       Number of pixels in Y direction
System.useGPU = 0;          % 1 or 0    Use GPU to accelerate computation. Effective when Nx, Ny is large (e.g. 600*800).
System.maxiter = 50;        % int       Number of iterations (for all methods explored)    
System.GSoffset = 1e-6;     % float>0   Regularization constant to allow low light background in 3D Gerchberg Saxton algorithms
System.source = ones(System.Nx, System.Ny); % amplitude of beam??
x1 = -System.Nx/2:System.Nx/2-1;   
y1 = -System.Ny/2:System.Ny/2-1;
RR = floor(System.Nx/2);
[x,y]=meshgrid(x1,y1);
circle = (x).^2+(y).^2;  
System.source(find(circle>RR^2))=0;  % diaphragm
System.focal_SLM = 0.4;   % unit:m   focal length of the telescope lens after slm.
System.focal_tube = 0.5;
System.focal_obj = 0.0072; % unit:m
System.ObjNA = 1.05;
System.ObjD = 15.12; %mm
System.ObjRI = 1.33;
Fx = System.lambda * System.focal_obj * System.focal_SLM / System.focal_tube / System.psSLM; % FOV, m
System.delta_x = Fx/System.Nx;

%% Target pattern definition
Depths = [0 0 0]*1e-6; % unit:m
System.tilt_posi_x = [175e-6,-330e-6,-50e-6];
System.tilt_posi_y = [-120e-6,-20e-6,330e-6];
tilt_x = round(System.tilt_posi_x/System.delta_x)+round(System.Nx/2);
tilt_y = round(System.tilt_posi_y/System.delta_x)+round(System.Nx/2);

% 1. multiple foci
I_target_focus = zeros(System.Nx,System.Nx,length(Depths));
for k = 1:length(tilt_x)
    I_target_focus(tilt_x(k),tilt_y(k),k)=1;
end

% 2. Beaded-ring pattern
r = 5e-6;            % radius of target pattern(s) (m) 
System.cor_x = x1*System.delta_x;
System.cor_y = y1*System.delta_x;
N=8;  % foci per beaded-ring pattern
I_target_BR = zeros(System.Nx,System.Nx,length(Depths));
for k = 1:length(tilt_x)
    for nn = 1:N
            xx =r*cos(2*nn*pi/N);
            yy =r*sin(2*nn*pi/N);
            particle_x =xx+System.tilt_posi_x(k);
            particle_y =yy+System.tilt_posi_y(k);
            [~,Index_x] = min(abs(System.cor_x-particle_x));
            [~,Index_y] = min(abs(System.cor_y-particle_y));
            I_target_BR(Index_x,Index_y,k) = 1;
    end
end

%% hologram generation
target_tem_focus = ones(length(Depths),1);
% target_tem_focus = [2 1 0.5]; % weight modified by real intensity
Iw_focus= sum(target_tem_focus(:))./(target_tem_focus+System.GSoffset)/length(Depths);
target_tem_BR = ones(length(Depths),1);
% target_tem_BR = [2 1 0.5]; % weight modified by real intensity
Iw_BR= sum(target_tem_BR(:))./(target_tem_BR+System.GSoffset)/length(Depths);
[HStacks] = function_Hstacks_cos(System,Depths);
[phase_focus]=WGS_para_noab_mod(System,I_target_focus,HStacks,Iw_focus);
[phase_BR]=WGS_para_noab_mod(System,I_target_BR,HStacks,Iw_BR);
%% pattern generation with SLM and save of hologram
sca
PsychDefaultSetup(2);
screens = Screen('Screens');
Screen('Preference','SkipSyncTests', 1);
screenNumber = 2; %SLM screen number
black = BlackIndex(screenNumber);
Screen('Preference', 'ConserveVRAM', 64);
Screen('Preference', 'SkipSyncTests', 1); % some internal self-tests and calibrations were skipped
[window, windowRect] = PsychImaging('OpenWindow',screenNumber,black); 
topPriorityLevel = MaxPriority(window);
ifi = Screen('GetFlipInterval', window)*2;
waitframes = 2;
Priority(topPriorityLevel);
vbl = Screen('Flip', window);
Nx=600; % pixel number of SLM
Ny=792; 
phase_tem = zeros(Nx,Ny);
% BR test
phase_tem(:,Ny/2-Nx/2:Ny/2+Nx/2-1) = phase_BR+pi;
phase_tem(Nx/2-System.Nx/2+1:Nx/2+System.Nx/2,Ny/2-System.Ny/2:Ny/2+System.Ny/2-1) = phase_BR+pi;
corrected_pat = double(imread('CAL_LSH0701883_1040nm.bmp')); % corrected pattern
wave_ratio = 211/255; % corrected coefficient of 1040 nm wavelength
phase_00 = double(uint8((phase_tem)/ 2 / pi * 255)) + corrected_pat;
phase_01 = mod(phase_00, 256)*wave_ratio;
phase_02 = uint8(phase_01);
imwrite(phase_02,'BR_ori.bmp','bmp')
Screen('PutImage', window, phase_02);
Screen('Flip', window);
save('ori_3pattern.mat')
% focus test
phase_tem(Nx/2-System.Nx/2+1:Nx/2+System.Nx/2,Ny/2-System.Ny/2:Ny/2+System.Ny/2-1) = phase_focus+pi;
corrected_pat = double(imread('CAL_LSH0701883_1040nm.bmp')); 
wave_ratio = 211/255; 
phase_00 = double(uint8((phase_tem)/ 2 / pi * 255)) + corrected_pat;
phase_01 = mod(phase_00, 256)*wave_ratio;
phase_02 = uint8(phase_01);
imwrite(phase_02,'focus_ori.bmp','bmp')
Screen('PutImage', window, phase_02);
Screen('Flip', window);





