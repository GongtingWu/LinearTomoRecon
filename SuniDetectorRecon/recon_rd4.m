% Revise on Feb.29 2016 again because Dr.Lu needs new images for his grant 
% proposal
%
% Make several changes, including output format, data truncation artifact
% correction. 


% This is the modified version of recon_d.m. It corrects the geometry
% misalignment
% 
% Gongting Wu, Oct.15th 2014 on UNC-CH


clear,
close all,
clc,

%--------------------------------------------------------------------------
%%    Source & Detector & Recon Space Geometry Input
%--------------------------------------------------------------------------
% Input directory
indir='D:\Personal Folder\Gongting Wu\____#Data#____\2015\0215 RMI dental phantom';

% Output directory
otdir='D:\Personal Folder\Gongting Wu\___#Result#___\2016\0229 RMI synthetic image_2';

% Actual height of the recon space
rh=24;

% Naming sequence
nseq=[7,8,10,12,14,15,17,18,21,22,24,26,27,29,30];

% Half projection mode
% hproj=[1,3,6,8,9,11,13,15];

% Original crop 
dl = 3;
dr = 3;
dt = 3;
db = 3;

% Pixel to crop
ct=2;
cb=2;
cl=2;
cr=2;

% pre-precessing factor
ppf=32768;

% Recon slice thickness (in cm)
thk=0.5;

% Constant used to amplified output image from float to uint32
otamp=10000;

% Source y distance (from N14(negative) to P14(positive))
% sy=[152.22,136.23,128.23,112.23,96.24,88.24,72.24,56.25,48.25,32.25,16.26,8.26,-7.74,...
%     -15.73,-31.73];
cd(indir)
load Suni-Geo.mat

Y_new=Y_new(nseq);
X_new=X_new(nseq);

% Number of tomo projections
% gnum=length(sy);
if(exist('hproj','var'))
    gnum=length(hproj);
else
    gnum=length(Y_new);
end


% Source to detector(image) distance in mm
sid=mean(Z_new);

% Object to detector distance
od=4;

% pixel size in x direction (row)
dxp=782;

% pixel size in y direction (column)
dyp=1066;

% pixel size in x direction (in mm)
dxs=0.033;

% pixel size in y direction (in mm)
dys=0.033;

% Crop the image
dxp=dxp-cl-cr;
dyp=dyp-ct-cb;

% detector length in x direction
dxl=dxp*dxs;

% detector length in y direction
dyl=dyp*dys;


% Height sampling of the recon space
rz=rh/thk;

% Correct for geometry misalignment
% Find the tilted angle of the linear source
sgeo=fit(X_new,Y_new,'poly1');
sang=acotd(-sgeo.p1);

% Correct for source_X and source_Y
sx = -sgeo.p2/sgeo.p1*cosd(sang);
sy = Y_new./cosd(sang);

% Half projection mode
if(exist('hproj','var'))
    nseq(setdiff(1:15, hproj))=[];
    sy(setdiff(1:15, hproj))=[];
end

%--------------------------------------------------------------------------
%%    Recon Parameters
%--------------------------------------------------------------------------
% Set regulation terms in the SART
opt.nonneg=1;
opt.box=1;

% Converge speed
opt.lambda='psi2mod';

% Number of iterations
itn=20;

% View direction: if vd=1, the reocon image will start from the bottom to
% the top, like recon in RTT; if vd=0, the recon image will start from top
% to the bottom.
vd=1;

% Set the system info. Actual height per z-pixel will be computed for each 
% recon slice
sysinfo.sh = sy; % Source y position (For each fan-beam 2D reconstruction)
sysinfo.sx = sx; % Source x position
sysinfo.sv = sid; % Source y position (vertical)
sysinfo.hsamp = dyp; % Horzontal sampling
sysinfo.vsamp = rz; % Vertical sampling
sysinfo.odd = od;% Object to detector distance
sysinfo.dpl = dyl;% Detector length
sysinfo.dpw = dxl; % Detector width
sysinfo.dp = dyp;% Number of pixels in detector
sysinfo.sdp = dxp;% Number of pixels in detector
sysinfo.rh = rh;% Actual height of the recon space

% Set the geo parameter to sample back the recon image
sxp=(sx+dxl/2)/dxs;
geo=single([sid,rh,od,sxp]);


%--------------------------------------------------------------------------
%%    Loads images, pre-processing images & Re-group data
%--------------------------------------------------------------------------

% Define the number of input dicom image, it should be named from 1 to N. 
disp('Input dicom image should indexed as 1,2,3,...10,11,12,...')

A=zeros(dyp,dxp,gnum);
for i=1:gnum
    tmp=dicomread([indir '\Image_' num2str(nseq(i),'%02d') '.dcm']);
    tmp=single(tmp)./ppf;       
    tmp(tmp>0.97)=1;
    
    % Correct the boundary
    tmp(1:dt,:)=repmat(tmp(dt+1,:),dt,1);
    tmp(end-db+1:end,:)=repmat(tmp(end-db,:),db,1);
    tmp(:,1:dl)=repmat(tmp(:,dl+1),1,dl);
    tmp(:,end-dr+1:end)=repmat(tmp(:,end-dr),1,dr);
    
    % Rotate, pad and crop the image
    tmp1 = padarray(tmp,[50,50],'replicate');
    
    tmp2=imrotate(tmp1,sang,'bicubic','crop');
    
    tmp3=tmp2(ct+51:end-cb-50,cl+51:end-cr-50);
    
    tmp3(tmp3<0.004)=0.004;
    
    A(:,:,i)=log(tmp3)*(-1);
end

% Regroup all data into fan beam projection
A=permute(A,[1,3,2]);
B=zeros(dyp*gnum,dxp);
for i=1:dxp
    tmp=A(:,:,i);
    B(:,i)=tmp(:);
end


% Release memory
clear A tmp

disp('Data is re-grouped into fan beam projection')


%--------------------------------------------------------------------------
%%    Reconstruction
%--------------------------------------------------------------------------
% Enable parallel computing
if isempty(gcp('nocreate'))
    parpool (feature('numCores'))
end

% Create recon space matrix
rb=zeros(rz,dyp,dxp);

nB=zeros(dyp*gnum,dxp);

% % Compute the parameter required to find system matrix A
% theta = atan(abs((1-1/2)*dxs-sx-dxl/2)./sid);
% 
% % Actual vertical length in the fan beam recon slice (z direction)
% hpz0 = 1/cos(theta);

% Compute the system matrix A
% A=fanbeam1(sysinfo,hpz0);
[Ax,Ay,wt] = ddProjector1(sysinfo);

tic;
parfor i=1:dxp
    % Compute the parameter required to find system matrix A
    theta = atan(abs((i-1/2)*dxs-sx-dxl/2)./sid);
    
    
    % Actual vertical length in the fan beam recon slice (z direction)
    hpz = 1/cos(theta);

    newA = sqrt(Ax + Ay*hpz^2);
    
    projData = B(:,i).*wt;
    
    % Recon the fanbeam slice using the system matrix A and the re-grouped
    % projection data B(:,i)
    sl=sart2(newA,projData,itn,zeros(size(newA,2),1),opt);

    % Generate a synthetic projection image
    nB(:,i)=newA*sl;
    
    % Reshape 1d recon array to 2d image and save in i-th slice
    sl=reshape(sl,rz,dyp);
    
    % Add recon slice image to the 3d recon image matrix
    rb(:,:,i)=sl;

    % Disp status
    disp(['Recon ' num2str(i) 'th slice...'])
end

% Permute the matrix to [x,y,z] for re-scaling
ra=permute(rb,[3,2,1]);

% Perform rescaling
r=sampb(ra,geo);
%r=ra;

% Permute back to [y,x,z] for display
r=permute(r,[2,1,3]);

toc;

% Turn off warning
warning('OFF','images:initSize:adjustingMag');

% Create output folder if it is exist
if(exist(otdir,'dir')==0)
    mkdir(otdir);
    mkdir([otdir '\Synthetic']);
end

% Save the recon result
for i=1:size(r,3)    
    if(vd==1)
        dicomwrite(uint16(imrotate(otamp*r(:,:,end+1-i),-sang,'bicubic','crop')),[otdir '\' num2str(i) '.dcm']);
    elseif(vd==0)
        dicomwrite(uint16(imrotate(otamp*r(:,:,i),-sang,'bicubic','loose')),[otdir '\' num2str(i) '.dcm']);
    end
end

% Save the synthetic projection image
nA=zeros(dyp,gnum,dxp);

for i=1:dxp
    tmp=reshape(nB(:,i),dyp,gnum);
    nA(:,:,i)=tmp;
end
nA=permute(nA,[1,3,2]);

for i=1:gnum
    dicomwrite(uint16(ppf*(1-exp(-nA(:,:,i)))),[otdir '\Synthetic\' num2str(i) '.dcm']);
end


disp('.');
disp('Reconstructed images are saved, this computer is yours nows :) ')


