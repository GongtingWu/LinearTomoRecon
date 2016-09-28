% This script includes all parameters used for reconstruction
% Gongting Wu, Sep.27 2016

%% I/O settings
dirHome = 'D:\Personal Folder\Gongting Wu\';
% Input directory
indir=[dirHome '____#Data#____\2016\0927 Dental\Teeth Images\'];
% Output directory
otdir=[dirHome '___#Result#___\2016\0927 Dental\Teeth Images-2\'];
% Constant used to amplified output image from float to uint16
otAmp=10000;

%% Pre-processing settings
% pre-precessing factor (2^16, 2^15, 2^14...)
ppf=9400;%2^14; %32768;
% Raw data thresholding
rawUp = 0.98;
rawLo = 0.0035;
% Raw data boundary correction 
dl = 3;
dr = 3;
dt = 3;
db = 3;

%% Reconstruction settings
% Recon space height (mm)
rh=24;
% Recon slice thickness (mm)
thk=0.5;
% Pixel to crop (Top, bottom, left, right)
% Flip for dental image from Corner
pixCrop = [160,10,10,10];
% Set regulation terms in the SART
opt.nonneg=1;
opt.box=1;
% Converge speed
reconLambda='psi2mod';
% Number of iterations
itn=20;
% Generate synthetic projection image
isSynImg = 0;

%% Geometry settings
disp(['Please make sure that order of the raw images matches that of ',...
    'the system geometry result']);

% Indicate the location of the geometry calibration file
geoFile = [dirHome, '\____#Data#____\2016\0927 Dental\geo.mat'];
load(geoFile);
geometryProfile = 'dental';

switch geometryProfile
    case 'dental'
        % Object to detector distance
        od=4;
        % pixel size in x,y direction (row)
        dxp=1440; % rows
        dyp=1920; % columns

        % pixel size in x, y direction (mm)
        dxs=0.0185;
        dys=0.0185;
        % Source positions
        sid=mean(z_locs);
        Y_new=x_locs - dxp*dxs/2;
        X_new=y_locs;
    case 'dbt'
        % TODO
        
        % Object to detector distance
        od=25;
    case 'dct'
        % TODO
end
