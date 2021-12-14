%% Load hologram (human spermozoids)
clear
close all
clc

load('holo_example')

%% DarkFocus
% Parameters
opts.dist = 0.085; % (mm)
opts.distRange = [-50,50]; % (um)
opts.distStep = 1; % (um)
opts.lambda = 0.45; % (um)
opts.pixSize = 5.5; % (um)
opts.mag = 18;
adv = [];

% Optional parameters
adv.ROI = 0;
adv.showProgress = 1;
adv.otherMethods = 0;

% Algorithm
[z,fC,DarkVolume] = DarkFocus(holo,opts,adv);
%% Simple GUI to show the results
DarkFocusGUI_single(fC, opts, DarkVolume)

%% Reconstruct darkfield hologram at z found by DarkFocus_v03
obj = AS_propagate(holo-imgaussfilt(holo,30),z*1000,1,opts.lambda,...
    opts.pixSize/opts.mag);
figure; imagesc(abs(obj)); colormap gray; axis image; title('amplitude')

%% Reconstruct hologram at z_new selected in DarkFocusGUI_single
% e.g., try z_new = 7 - a single spermozoid is in focus, while all others
% are out of focus
obj = AS_propagate(holo-imgaussfilt(holo,30),opts.dist*1000+z_new,1,opts.lambda,...
    opts.pixSize/opts.mag);
figure; imagesc(abs(obj)); colormap gray; axis image; title('amplitude')

%% Simple GUI to show the results at 3 different planes
DarkFocusGUI_triple(fC, opts, DarkVolume)

%% Reconstruct hologram at outDF selected in DarkFocusGUI_triple
obj1 = AS_propagate(holo-imgaussfilt(holo,30),opts.dist*1000+outDF(1),1,opts.lambda,...
    opts.pixSize/opts.mag);
figure; imagesc(abs(obj1)); colormap gray; axis image; title('amplitude 1')
obj2 = AS_propagate(holo-imgaussfilt(holo,30),opts.dist*1000+outDF(2),1,opts.lambda,...
    opts.pixSize/opts.mag);
figure; imagesc(abs(obj1)); colormap gray; axis image; title('amplitude 2')
obj3 = AS_propagate(holo-imgaussfilt(holo,30),opts.dist*1000+outDF(3),1,opts.lambda,...
    opts.pixSize/opts.mag);
figure; imagesc(abs(obj3)); colormap gray; axis image; title('amplitude 3')

%% Auxiliary function

function uout = AS_propagate(uin, z, n0, lambda, dx)

[Ny,Nx] = size(uin);
dfx = 1/Nx/dx;
dfy = 1/Ny/dx;
fx=(-Nx/2:Nx/2-1)*dfx;
fy=(-Ny/2:Ny/2-1)*dfy;
[FX,FY] = meshgrid(fx,fy);

%FT of the input field
FTu = fftshift(fft2(uin));

% %generation of the transfer function
FZ = sqrt((n0/lambda).^2-FX.^2-FY.^2);
FZ(~isreal(FZ))=0;
TF = exp(1i*2*pi*z*FZ);

% multiplication with the transfer function
FTu = FTu.*TF;

% inverse FT
uout = ifft2(fftshift(FTu));
end