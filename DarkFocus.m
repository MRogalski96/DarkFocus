function [z,fC,DarkVolume] = DarkFocus(holo,opts,adv)
% Autofocusing algorithm for lensless digital in-line holographic 
% microscopy (DIHM). It finds the distance between the hologram (holo) and 
% object focus plane (z, [mm]). It also returns the DarkFocus focus curve 
% (fC) and the reconstructed DarkVolume (result of darkfield hologram 
% propagation to investigated planes)
%
% Inputs:
%   holo - lensless microscopy hologram
%   opts - reconstruction options
%       opts.dist - propagation distance (mm) - a distance between camera
%           and the middle of object plane (initial guess of the z)
%       opts.distRange - tested range of the propagation distances (um)
%           opts.distRange = [lower bound, upper bound]; (see drawing
%                            below)
%       opts.distStep - sampling in z direction (um)
%           propagationRange = opts.dist*1000+opts.distRange(1) : ...
%               opts.distStep : opts.dist*1000+opts.distRange(2)
%       opts.lambda - wavelength (um)
%       opts.pixSize - camera pixel size (um)
%       opts.mag - magnification of the microscope setup
%       opts.n0 - background refractive index (default = 1)
%   adv - advanced options (optional)
%       adv.ROI - region of interest that will be investigated,
%           adv.ROI = [x0, y0, x_length, y_length];
%               e.g., adv.ROI = [1,1,500,500];
%           if adv.ROI = 1; - you will choose the ROI on the reconstructed
%               holo(:,:,1) at opts.dist distance (mark the ROI and dobule
%               click inside it)
%           if adv.ROI = 0; - full image will be investigated (default)
%       adv.showProgress - display algorithm progress (every 10%) in
%                          command wintow
%           adv.showProgress = 0; - no 
%           adv.showProgress = 1; - yes (default)
%       adv.backRemov - background removing method
%           adv.backRemov = 1; - hologram background is obtained through
%               gaussian filtering
%           adv.backRemov = 2D matrix - adv.backRemov is set as the 
%               hologram background 
%       adv.otherMethods - calculate the focus curves for other methods
%           adv.otherMethods = 0 - no, only DarkFocus (fC(1,:)) (default)
%           adv.otherMethods = 1 - yes, calculate also Tamura (fC(2,:)) 
%                                  and Dubois (fC(3,:)) methods
%
% Outputs:
%   z - distance from the focus plane to the camrea plane (mm)
%   fC - DarkFocus focus curve (also Tamura and Dubois focus curves if 
%        adv.otherMethods == 1)
%   DarkVolume - a set of generated darkfield amplitudes
%
% Ilustrative drawing:
% dR(2) dR(1)   - here dR = opts.distRange; dR(1) < 0 and dR(2) > 0
% <---|<--|
%  _______
% |   |   |              |
% |tested |              |
% |prop.  |              |detector
% |range  |              |plane
% |___|___|              |
%      <-----------------|
%            opts.dist
%
% Created by:
%   Mikołaj Rogalski,
%   mikolaj.rogalski.dokt@pw.edu.pl
%   Institute of Micromechanics and Photonics,
%   Warsaw University of Technology, 02-525 Warsaw, Poland
%
% Last modified: 14.12.2021

%% Deal with the input
if nargin < 2
    error('Not enaugh input arguments ("holo" and/or "opts")');
end
if ~isfield(opts,'n0'); opts.n0 = 1; end
if nargin < 3; adv = []; end
if ~isfield(adv,'ROI'); adv.ROI = 0; end
if ~isfield(adv,'showProgress'); adv.showProgress = 1; end
if ~isfield(adv,'backRemov'); adv.backRemov = 1; end
if ~isfield(adv,'otherMethods'); adv.otherMethods = 0; end
%% Initial processing
% Propagation range
propRang = opts.dist*1000+opts.distRange(1):opts.distStep:...
    opts.dist*1000+opts.distRange(2);
% Pixel size in object plane
dPix = opts.pixSize/opts.mag;
% Processed hologram
holo = double(holo);
% Calculating hologram background
if adv.backRemov == 1
    bckr = imgaussfilt(holo,30);
else
    bckr = adv.backRemov;
end
holo_dark = holo - bckr;

% ROI selection
ROI = adv.ROI;
if length(ROI) == 1
    if ROI == 1
        Obj = AS_propagate(holo_dark,opts.dist*1000,1,opts.lambda,dPix);
        ROI0 = ROISelection(abs(Obj));
%         ROI0 = ROISelection(holo);
        holo_dark = holo_dark(ROI0(2):ROI0(2)+ROI0(4)-1,...
            ROI0(1):ROI0(1)+ROI0(3)-1);
        holo = holo(ROI0(2):ROI0(2)+ROI0(4)-1,ROI0(1):ROI0(1)+ROI0(3)-1);
    end
elseif length(ROI) == 4
    holo_dark = holo_dark(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
    holo = holo(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
end

[Sy,Sx] = size(holo_dark);
Sz = length(propRang);

%% Focus curve(s) calculation
DarkVolume = zeros(Sy,Sx,Sz);

DF = zeros(1,length(propRang));
pB = 0;
if adv.otherMethods == 1
    ToG = zeros(1,length(propRang));
    Dubois = zeros(1,length(propRang));
end

for mm = 1:Sz
        Obj1 = AS_propagate(holo_dark,propRang(mm),1,opts.lambda,dPix);
        % DarkFocus
        amp = abs(Obj1);
        DarkVolume(:,:,mm) = amp;
        [gx, gy] = gradient(amp);
        gradFilt = gx.^2+gy.^2;
        DF(mm) = var(gradFilt(:));
    
    if adv.otherMethods == 1
        Obj0 = AS_propagate(holo,propRang(mm),1,opts.lambda,dPix);

        % Tamura
        amp0 = abs(Obj0);
        [gx,gy] = gradient(amp0);
        grad=gx.^2+gy.^2;
        ToG(mm)=sqrt(sqrt(var(grad(:))))./sqrt(mean2(grad));
        
        % Dubois
        rre = real(Obj0) - imgaussfilt(real(Obj0),10);
        iim = imag(Obj0) - imgaussfilt(imag(Obj0),10);
        Dubois(mm) = sum(sum(abs(rre + 1i*iim)));
    end
    if adv.showProgress == 1
        p = floor(mm/Sz*10);
        if p ~= pB
            disp(['DarkFocus; Completed ',num2str(round(mm/Sz*100,1)),'%'])
        end
        pB = p;
    end
end

% normalization 0-1
DF = DF - min(DF); DF = DF/max(DF);
fC(1,:) = DF;
[~,loc] = max(fC);
z = propRang(loc)/1000;

if adv.otherMethods == 1
    ToG = ToG - min(ToG); ToG = ToG/max(ToG);
    fC(2,:) = ToG;
    Dubois = Dubois - min(Dubois); Dubois = Dubois/max(Dubois);
    fC(3,:) = Dubois;
end

end

%% Auxiliary functions
function Position = ROISelection(in)

% % 'Position',[x,y,w,h] [left upper width height]

figure(100);
img = imagesc(in); colormap gray; 
title('Select ROI and then double-click inside');
[X,Y,I2,Position] = imcrop(img);
Position = round(Position);
if ( mod(Position(3),2) == 1 ), Position(3) = Position(3)-1; end
if ( mod(Position(4),2) == 1 ), Position(4) = Position(4)-1; end
% out = in(Position(2):Position(2)+Position(4)-1,Position(1):Position(1)+Position(3)-1);
delete(gcf);
end

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