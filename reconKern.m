% Reconstruction loop

% Enable parallel computing
if isempty(gcp('nocreate'))
    parpool (feature('numCores'));
end

% Allocate recon space matrix
rb=zeros(rz,dyp,dxp,'single');
nB=zeros(dyp*gnum,dxp, 'single');

parfor i=1:dxp
    % Compute the parameter required to find system matrix A
    theta = atan(abs((i-1/2)*dxs-sx-dxl/2)./sid);
    
    % Actual vertical length in the fan beam recon slice (z direction)
    hpz = 1/cos(theta);

    newA = sqrt(Ax + Ay*hpz^2);
    
    % Truncation data correction
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