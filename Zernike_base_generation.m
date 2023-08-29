%%% Zernike base generation
%%% JC,2022/6/28

% SLM
Nx = 600;
Ny = 792;
% obj
refractive_n = 1.33;
NA = 1.05;
Diameter = 15; % mm,obj aperture
Diameter_act = 15; % mm,actual aperture

% Zernike
Zernike_num = 5; % Zernike mode 
n = zeros(Zernike_num,1);
m = zeros(Zernike_num,1);
for j=0:Zernike_num  
   n(j+1)=ceil((-3+sqrt(9+8*j))/2);   %highest power or order of the radial polynomial term 
   m(j+1)=2*j-n(j+1)*(n(j+1)+2);                %azimuthal frequency of the sinusoidal component 
end  
del_num = [1,2,3,5]; %delete constant, tip, tilt and defocus
n(del_num) = [];
m(del_num) = [];
X = linspace(-Diameter_act/2,Diameter_act/2,Nx);
Y = linspace(-Diameter_act/2,Diameter_act/2,Nx);

phase_Zer = zeros(Nx,Nx,length(n));
for j = 1:length(n)
    n_tem = n(j);
    m_tem = m(j);
    phase_Zer(:,:,j) = zernike_mod(n_tem,m_tem,X,Y,Diameter,Diameter_act);
end
save('base_Zernike_phase_upto27_fullpupil_withouttilt..mat','phase_Zer','n','m');
