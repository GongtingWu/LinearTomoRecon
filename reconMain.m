%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               AFVR One v1.0
%
% Author: Gongting Wu  |  Zhou Lab, UNC Chapel Hill

% 

% Clear all info
clear,
close all,
clc,

% Load recon parameters and process the data
reconProc;

% Compute the system matrix A
% A=fanbeam1(sysinfo,hpz0);
[Ax,Ay,wt] = ddProjector(sysinfo);

tic;
% Reconstruction Loop
reconKern;
toc;

% Permute the matrix to [x,y,z] for re-scaling
ra=permute(rb,[3,2,1]);

% Perform rescaling
r=sampb(ra,geo);
% r=ra;

% Permute back to [y,x,z] for display
r=permute(r,[2,1,3]);
toc;

%% Save reconstruction images
% Turn off warning
warning('OFF','images:initSize:adjustingMag');
% Save the recon result

if(strcmp(geometryProfile, 'dental'))
    for i=1:size(r,3)    
        dicomwrite(uint16(imrotate(otAmp*r(:,:,end+1-i),-sang,...
            'bicubic','crop'))',[otdir '\Image_' num2str(i,'%02d') '.dcm']);
    end
else
    for i=1:size(r,3)    
        dicomwrite(uint16(imrotate(otAmp*r(:,:,end+1-i),-sang,...
            'bicubic','crop')),[otdir '\Image_' num2str(i,'%02d') '.dcm']);
    end
end

if(isSynImg)
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
end

disp('Reconstruction completed!')


