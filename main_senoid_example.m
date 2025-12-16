

% Author: Vinicius Soares de Angelis   
% email: vinicius.angelis@usp.br 

% Conformal holography with curved light sheets
% https://doi.org/10.1364/OE.536859


close all;
clear all;
clc;
%% input parameters

% target (unfolded) profile

Ao1 = imread('template_shapes.png');
Ao1 = imrotate(Ao1,90);

% input parameters
                            
M = 121; % number of linear FWs ( = number of pixels in x direction desired for the the input image)
Mz = 201; % number of pixels in z direction desired for the the input image
          
lambda0 = 532e-9; %[m] operating wavelength of the RGB laser in free space


nm = 1; % refractive index of the medium (real positive value)
r0 = 35e-6; % [m] spot size radius of each FW
L_im = 500e-3; % [m] longitudinal distance of each FW ( = image length in z direction)
u = 0; % order of the Bessel beams


fprintf('Number of linear FWs: %d \n',M);
fprintf('Spot size radius of each FW: %f um\n',r0*1e6);
fprintf('Longitudinal distance of each FW: %f mm\n',L_im*1e3);

Ao = imresize(Ao1, [M Mz]); % resize input image
A = rgb2gray(Ao);
imshow(A)

syms kz z rho phi
digits(100);

I1 = single(A); % converts gray scale values to float
I1 = I1/max(max(I1)); % normalization of input image intensity 
u1 = sqrt(I1); % normalized field amplitude of input image  
figure
imshow(u1)
   
w0 = (2*pi*Const.c0)/lambda0; %[rad/s] operating angular frequency
lambda = lambda0/nm; %[m] wavelength of the medium
k = nm*(2*pi/lambda0); %[rad/m] wavenumber of the medium

 
L = 2*L_im; % the field must be designed only in 0 < z < L/2 for a non-absorbing medium  

Q = k*(1-0.5*(2.4048/(r0*k))^2);
Nmax = floor((L/(2*pi))*(k-Q));
krho0 = sqrt(2)*sqrt(1-Q/k)*k;
q = k/(L*krho0);
N = Nmax;

kzc = Q-2*pi*N/L; % longitudinal wave number of the Bessel beam with the highest axicon angle
axicon_max = acos(kzc/k); % highest axicon angle
krho_max = sqrt(2)*sqrt(1-(kzc/k))*k; % transverse wave number of the Bessel beam with the highest axicon angle  
Rab = L_im*sqrt((k/kzc)^2-1); % minimum radius aperture to generate each FW

fprintf('Guassian beam waist (1/q): %f mm \n',(1/q)*1e3);
fprintf('Number of Bessel beams in each FW superposition (2 N+1): %d \n',2*N+1);
fprintf('Paraxiallity level (Q/k): %f \n',Q/k);
fprintf('Highest axicon angle (º): %f \n',axicon_max*(180/pi));
fprintf('Corresponding longitudinal wavenumber (kzc/k): %f \n',kzc/k);
fprintf('Minimum pixel size (um): %f \n',0.5*1e6/krho_max);
fprintf('Minimum radius aperture each FW (mm): %f \n',Rab*1e3);

%% design of the curved surface
% here a sinonoidal curve (Fig. 2c)

Xmax = 1.5e-3;
Ymax = 0.25*Xmax;
freq = 1;
x = linspace(-Xmax,Xmax,301);
y = Ymax*sin(freq*pi*x/Xmax);

p = [x',y'];
q0 = curvspace(p,M);

x0 = q0(:,1);
y0 = q0(:,2);


figure
plot(x0*1e3,y0*1e3,'ko')
xlabel('x (mm)')
ylabel('y (mm)')


d = sqrt((x0(1:end-1)-x0(2:end)).^2+(y0(1:end-1)-y0(2:end)).^2);
% separation distances betwwen light threads

Zmax = L_im;
Zmin = 0;
Quant_z = Mz;
zp = linspace(Zmin,Zmax,Quant_z);
step = zp(2)-zp(1);


Lxab = 2*Rab+2*Xmax; % minimum dimensions of a rectangular aperture (placed at z = 0) able to generate
Lyab = 2*Rab+2*Ymax;                 % the set of FWs   

fprintf('Minimum rectangular aperture dimension in x (mm): %f \n',Lxab*1e3);
fprintf('Minimum rectangular aperture dimension in y (mm): %f \n',Lyab*1e3);

fprintf('percentage of SLM dim in x : %f \n',100*Lxab/15.4e-3);
fprintf('percentage of SLM dim in y : %f \n',100*Lyab/9.6e-3);

%% compute the superposition of the light threads
% field stored in a matlabFunction Ex (x-polarized light)
Ex = 0;
tic
for i=1:length(x0)
    aux = sum(u1(i,:).*function_gate(zp-step/2,zp+step/2));
    %aux = function_gate(Zmin-10,Zmax+10);
    [An,in,n] = profile_fwbg(aux,L,N);
    rho0 = sqrt(x0(i).^2+y0(i).^2);
    phi0 = atan2(y0(i),x0(i));
    aux_x = field_fwbg(u,1,An,in,n,L,k,Q,q,rho0,phi0);
    Ex = Ex + aux_x;
    
    fprintf('x0m(%d): %f mm \n',i,x0(i)*1e3);   
end

fEx = matlabFunction(Ex,'vars', [rho phi z]);
toc

%% Evaluate the field at the curved surface  
tic
Qz = 201;
xx = interpn(x0,2);
zz = linspace(Zmin,Zmax,Qz);
[XX,ZZ] = meshgrid(xx,zz);
YY = Ymax*sin(freq*pi*XX/Xmax);
RHO = sqrt(XX.^2 + YY.^2);
PHI = atan2(YY,XX);


E2E0 = abs(fEx(RHO,PHI,ZZ)).^2;

rgbValues = spectrumRGB(lambda*1e9);
discr_color = 255;
cmap = [linspace(0,1,discr_color)'*rgbValues(1) linspace(0,1,discr_color)'*rgbValues(2) linspace(0,1,discr_color)'*rgbValues(3)];

figure
surf(ZZ*1e3,-XX*1e3,-YY*1e3,(E2E0));grid
shading flat
colormap(cmap)
colorbar 
h = gca;
h.Visible = 'On';
axis square
shading flat
xlabel('z (mm)')
ylabel('x (mm)')
zlabel('y (mm)')
set(gca,'FontSize',15,'FontWeight','bold')
toc

%% Encoding of the field at z = 0 into a CGH 
% 
tic
n_x = 1920;  %Number of pixels in the horizontal direction x
n_y = 1200;  %Number of pixels in the vertical direction y
SS = 1 ;
xmax =  0.0077  ; %Number of pixels multiplied by the pixel pitch (8 um) divided by 2
xmin =  -xmax   ;
ymax =  0.0048  ; %Number of pixels multiplied by the pixel pitch (8 um) divided by 2
ymin =  -ymax   ;
nx = 960;
ny = 600;
xxa    =  linspace(xmin,xmax,nx); % The field is undersampled by a factor of 2 (we use 960 steps instead of 1920)
yya    =  linspace(ymin,ymax,ny); % The field is undersampled by a factor of 2 (we use 600 steps instead of 1200)
% Undersampling the field is equivalent to doubling the transverse spatial frequency
% krho, this helps generate the output pattern over a shorter propagation distance, making the measurement easier. 

[XXa,YYa] = meshgrid(xxa,yya);
rhorho=sqrt(XXa.^2+YYa.^2);
phiphi = atan2(YYa,XXa);
Zp = 0;
Exap = fEx(rhorho,phiphi,Zp);

[Det] =  GenerateHologram(Exap,nx,ny,XXa/2,YYa/2);

% 'GenerateHologram' function implements the CGH of type 3 
% of V. Arrizón (2007) - Ref. 38


% Normalize the phase function  
SLM1=Det/max(Det(:))*255;
H1 = SLM1;
% The aperture size of the SLM is 1920 X 1200
S1 = zeros(1200,1920);
% Assign the hologram to the center of the SLM
S1(301:900, 481:1440)    = H1;

% Digitize the hologram into 256 levels 
H = uint8(((round(S1*255./max(max(S1)))))) ;
%  
% Save the hologram into png file to be addressed onto SLM 
imwrite(H,('example_senoid.png'),'png') 

toc
%% Simulating the CGH 
% passing through the 4f sytem, selecting the first diff. order
% applying Fresnel propagation to get the profile at the 
% curved surface

figure
surf(XXa*1e3,YYa*1e3,abs(Exap).^2);grid
shading flat
colormap jet
colorbar 
xlabel('x (mm)')
ylabel('y (mm)')
title(['Squared amplitude    z = ' num2str(Zp*1e3),' mm']);
view([0 90])
xlim([xmin(1)*1e3 xmax(end)*1e3])
ylim([ymin(1)*1e3 ymax(end)*1e3])

F1 = ifftshift(fft2(fftshift(exp(1i*Det))));

dfx = 1/(2*nx*8e-6*0.5);
dfy = 1/(2*ny*8e-6*0.5);

fx = (-nx/2:nx/2-1)*dfx*2*pi;
fy = (-ny/2:ny/2-1)*dfy*2*pi;

[Fx,Fy] = meshgrid(fx,fy);
figure
surf(Fx,Fy,abs(F1).^2);grid
shading flat
colormap jet
colorbar 
xlabel('k_x')
ylabel('k_y')
title('Fourier plane');
view([0 90])
set(gca, 'CLim', [0, 10000]);

gy=300/(2*nx*8e-6);
gx=300/(2*ny*8e-6);

gyp = (gx/dfx);
gxp = (gy/dfy);

xc = (0.5*nx+gyp);
yc = (0.5*ny+gxp);
rdmax = sqrt(gyp^2+gxp^2);
r_krhomax_x = ((2*krho_max/(2*pi))/dfx);
r_krhomax_y = ((2*krho_max/(2*pi))/dfy);
alpha = 1;
rdx = round(r_krhomax_x*alpha);
rdy = round(r_krhomax_y*alpha);

Ff = F1.*((sqrt((Fx-gx*2*pi).^2+(Fy-gy*2*pi).^2))<2*krho_max*alpha);
F2 = Ff(round(yc-rdy):round(yc+rdy),round(xc-rdx):round(xc+rdx));

F2b = padarray(F2,[floor((ny-size(F2,1))*0.5)+1,floor((nx-size(F2,2))*0.5)+1],0,'both');
F2a = imresize(F2b,[ny nx]);

figure
image(abs(F2a).^2)

F3 = ifftshift(ifft2(fftshift(F2a)));

figure
surf(0.5*XXa*1e3,0.5*YYa*1e3,abs(F3).^2);grid
shading flat
colormap jet
colorbar 
xlabel('x (mm)')
ylabel('y (mm)')
title(' z = 0');
view([0 90])

tic
zt = (0.25*L_im);
nz = 200;
z = linspace(0,zt,nz);
x01 = interpn(x0,2);
npx = size(x01,2);
y01 = interpn(y0,2);
Up = zeros(nz,npx);

figure
for jj = 1:length(z)

dx = max(xxa)/nx;
dy = max(yya)/ny;
[U,u,v]=Fresnel(F3,lambda,z(jj),dx,dy,0);


surf(0.5*XXa*1e3,0.5*YYa*1e3,abs(U).^2);grid
shading flat
colormap(cmap)
colorbar 
xlabel('x (mm)')
ylabel('y (mm)')
title(' z = 0');
view([0 90])
title(['Fresnel propagation   z = ' num2str(z(jj)*1e3),' mm']);
pause(0.1);

x0p = 0.5*nx+round(0.5*x01/dx);
y0p = 0.5*ny+round(0.5*y01/dy);

    for j=1:npx
        Up(jj,j) = U(y0p(j),x0p(j)); 
    end
end
toc
%%
% profile at the curved surface
% simulated from the CGH
% this is what is expected to be obtain in the lab

[Xp,ZZp1] = meshgrid(x01,z);
[Yp,ZZp] = meshgrid(y01,z);

figure
surf(ZZp*1e3,0.5*Xp*1e3,0.5*Yp*1e3,(abs(Up).^2));grid
shading flat
colormap(cmap)
colorbar 
h = gca;
h.Visible = 'On';
axis square
shading flat
xlabel('z (mm)')
ylabel('x (mm)')
zlabel('y (mm)')
set(gca,'FontSize',15,'FontWeight','bold')
%%




