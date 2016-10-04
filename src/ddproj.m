function [A, l] = ddproj(geom)
% Optimized version of distance driven projector used in AFVR. All x-ray 
% path lengths are output as a nx1 vector *lx*, which are used to compute
% the specific weighting matrix for each fan volume planes. 
%
% Note: Trunation artifact weighting is automatically applied in the system
% matrix *A*.
%
% Author: Gongting Wu, Zhou Lab @ UNC Chapel Hill
% Date: Oct.4th 2016


% Display version
version = 'dd_1.10';
disp(['This is the distance driven version: ' version]);

% Discretization intervals in each dimension
srcX=geom.sh; % Source x position (horizontal)
sid=geom.sv; % Source y position (vertical)
reconSlice=geom.vsamp; % Vertical sampling
objDetDistance=geom.odd;% Object to detector distance
detLength=geom.dpl;% Detector length
numDetPixel=geom.dp;% Number of pixels in detector
reconHeight=geom.rh;% Actual height of the recon space

% ------------------------------------------------------------
%  Calculating all parameters for the forward projection loop
% ------------------------------------------------------------

% nA denotes the number of projections.
numProj = length(srcX);

% Length for each pixel
pixelSize=detLength/numDetPixel;

% Heigh for each voxel
voxelHeight=reconHeight/reconSlice;

% Initialize vectors that contains the row numbers, the column numbers and
% the values for creating the matrix A effiecently.
% rows = zeros((nx+ny)*nA*Dp,1);
rows = zeros(2*numReconSlice*numProj*numDetPixel,1);
cols = rows;
vals = rows;
idx = 0;

% No need to weight, direct weight
% % Initialize the weighting for the measurement
% wt = ones (numProj*numDetPixel,1);

for i = 1:numProj
    
    % The coordinate of the source location for this projection
    src = srcX(i);

    % Source to detector edge distance
    src_dist = abs(src) - detLength/2;

    % Computing the weighting factor
    if src_dist>0
        % Edge of the recon space rectangular
        top_edge = objDetDistance+reconHeight-voxelHeight/2;
        bot_edge = objDetDistance+voxelHeight/2;

        % The calculated distance should divide by the unit length "pixelSize" to
        % find the pixel location
        det_end = floor(top_edge/(sid-top_edge)*src_dist/pixelSize)+1;
        det_srt = ceil(bot_edge/(sid-bot_edge)*src_dist/pixelSize)+1;
        ind = det_srt:det_end;
        det_dist = (ind-1/2)*pixelSize;
        edge_height = det_dist./(src_dist + det_dist)*sid;
        w2 = (edge_height - objDetDistance)/reconHeight;

        if src<0
            idx1 = (i-1)*numDetPixel;
            wt(idx1+1:idx1+det_srt-1)=0.01;
            wt(idx1+det_srt:idx1+det_end)=w2;

        elseif src>0
            idx1 = i*numDetPixel;
            wt(idx1-det_srt-1:idx1)=0.01;
            wt(idx1-det_end:idx1-det_srt)=fliplr(w2);
        end
    end

    
    %% Start
    if abs(src)<= (detLength/2)
        for j = 1:numReconSlice
            slicePos = (j-1/2)*voxelHeight+(sid-reconHeight-objDetDistance);
            ratio = slicePos/sid;
            startPos = (src + detLength/2) * (1-ratio);
            startPoint = startPos/pixelSize;
            for k = 1:numDetPixel
                nsp = startPoint+(k-1)*ratio;
                nep = nsp + ratio;
                dp = (i-1)*numDetPixel+k;
                if(ceil(nsp) == ceil(nep))
                    idx = idx +1;
                    rows(idx) = dp;
                    vals(idx) = 1;
                    cols(idx) = floor(nsp)*numReconSlice + j;
                else
                    idx = idx +1;
                    rows(idx) = dp;
                    vals(idx) = (ceil(nsp) - nsp)/ratio;
                    cols(idx) = floor(nsp)*numReconSlice + j;
                    idx = idx +1;
                    rows(idx) = dp;
                    vals(idx) = (nep - floor(nep))/ratio;
                    cols(idx) = floor(nep)*numReconSlice + j;
                end
            end
        end
    elseif src<= (-detLength/2)
        for j = 1:numReconSlice
        slicePos = (j-1/2)*voxelHeight+(sid-reconHeight-objDetDistance);
        ratio = slicePos/sid;
        % Determine the start point
        startDetPos = (src+detLength/2)*(1-1/ratio);
        startPixel = startDetPos/pixelSize;
        
        % Determine if the start point is right on the pixel boundary
        if(ceil(startPixel) == floor(startPixel))
            startPixel = ceil(startPixel)+1;
            startPoint = 0;
        else
            startPoint = (ceil(startPixel)-startPixel)*ratio;
            idx = idx +1;
            rows(idx) = (i-1)*numDetPixel+ceil(startPixel);
            vals(idx) = startPoint/ratio;
            cols(idx) = floor(startPixel)*numReconSlice + j;
            startPixel = ceil(startPixel)+1;     
        end
        
        % Loop over all pixels in the detector
        for k = startPixel:numDetPixel
            nsp = startPoint+(k-startPixel)*ratio;
            nep = nsp + ratio;
            dp = (i-1)*numDetPixel+k;
            if(ceil(nsp) == ceil(nep))
                idx = idx +1;
                rows(idx) = dp;
                vals(idx) = 1;
                cols(idx) = floor(nsp)*numReconSlice + j;
            else
                idx = idx +1;
                rows(idx) = dp;
                vals(idx) = (ceil(nsp) - nsp)/ratio;
                cols(idx) = floor(nsp)*numReconSlice + j;
                idx = idx +1;
                rows(idx) = dp;
                vals(idx) = (nep - floor(nep))/ratio;
                cols(idx) = floor(nep)*numReconSlice + j;
            end
        end
        end
    elseif src>= (detLength/2)
        for j = 1:numReconSlice
            slicePos = (j-1/2)*voxelHeight+(sid-reconHeight-objDetDistance);
            ratio = slicePos/sid;
            
            startPos = (src + detLength/2) * (1-ratio);
            startPoint = startPos/pixelSize;
            
            % Find the last pixel
            endPos = detLength - (src - detLength/2) * (1/ratio-1);
            endPixel = endPos/pixelSize;
            
            % Determine if the end point is right on the pixel boundary
            if(ceil(endPixel) == floor(endPixel))
                endPixel = floor(endPixel)-1;
            else
                endPoint = endPixel-floor(endPixel);
                idx = idx +1;
                rows(idx) = (i-1)*numDetPixel+ceil(endPixel);
                vals(idx) = endPoint;
                cols(idx) = floor(endPixel)*numReconSlice + j;
                endPixel = floor(endPixel)-1;     
            end
            
            % Loop over all pixels in the detector
            for k = 1:endPixel
                nsp = startPoint+(k-1)*ratio;
                nep = nsp + ratio;
                dp = (i-1)*numDetPixel+k;
                if(ceil(nsp) == ceil(nep))
                    idx = idx +1;
                    rows(idx) = dp;
                    vals(idx) = 1;
                    cols(idx) = floor(nsp)*numReconSlice + j;
                else
                    idx = idx +1;
                    rows(idx) = dp;
                    vals(idx) = (ceil(nsp) - nsp)/ratio;
                    cols(idx) = floor(nsp)*numReconSlice + j;
                    idx = idx +1;
                    rows(idx) = dp;
                    vals(idx) = (nep - floor(nep))/ratio;
                    cols(idx) = floor(nep)*numReconSlice + j;
                end
            end           
        end
    end
end

% Truncate excess zeros.
rows(idx+1:end) = [];
cols(idx+1:end) = [];
vals(idx+1:end) = [];
vals = vals/numReconSlice;

% Create sparse matrix A from the stored values.
A = sparse(rows,cols,vals,numDetPixel*numProj,numReconSlice*numDetPixel);