function [r]=sampb(rb,geo)
% sampb function take 3d recon image matrix, and samples it back to its
% true aspect ratio
% 
% --------------Update on Apr.30th 2015--------------
%       1. Fix "size does not match" error
%
% --------------Update on Sep.15th 2014--------------
%           Now it support negative value 
%       (Source at the left of the detector: sx<0)

% Extract data from input
sd=geo(1); % Source to detector distance
rh=geo(2); % actual height of the recon space
od=geo(3); % object to detector distance
sx=geo(4); % x position of the source with the unit of pixel length related to the edge of the recon space

% Find recon image sampling
[rx,ry,rz]=size(rb);

% Height per recon slice
ph=rh/rz;

% Divide the problem into two case
if sx>rx
    error(['Sample back code only works when the source X position',...
        ' is within the detector. Need to modify the Sample back code!']);
end


% This is actually the true pixel size at each recon slice
plz=1-(od+((rz:-1:1)-1/2).*ph)./sd;

% Compute the starting point for each x-y plane slice
sp = (1-plz).*sx;% sx here is source x position plus half detector  
% (04/30/15) sx is the number of pixels from the left edge of the detector

% Compute the ending points
ep = (rx-sx)*plz+sx;

% First index in normal space 
fsp = round(sp); % Last point where no value is assigned due to the limited angle geometry
% Last index in normal space
fep = round(ep);
% Number of each interpolation points
itp_num = fep - fsp;

% First pixel position at each reconstruction slices
dist = (fsp+1-sp)./plz;

% Initialize r matrix
r=zeros(size(rb));
% Pad zeros
for i=1:rz
    r(1:fsp(i),:,i)=0;
    r(fep(i)+1:end,:,i)=0;
end
% 
% % 
% % 
% % % First index in the re-sampled space
% % fpr=(ceil(sp)-sp)./plz+1; 
% % % Last point of the normal space in re-sampled space
% % lpr = (fep-2-sp)./plz+1;

% Interpolate the value (Scale down to the correct geometry)
for j=1:ry
    for i=1:rz
        r(fsp(i)+1:fep(i),j,i)=interp1(1:rx,rb(:,j,i),((1:itp_num(i))./plz(i)+dist(i))');
    end
end