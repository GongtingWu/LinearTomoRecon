function [Ax,Ay,wt] = ddProjector1(info)
% Distance driven projector for the flat panel tomosynthesis reconstruction

% Assume the pixel size in each reconstruction slice is same as the
% detector pixel size

% Gongting Wu, Oct.3rd, 2015

version = 'dd_1.05';

disp(['This is the distance driven version: ' version]);

% Discretization intervals in each dimension
srcX=info.sh; % Source x position (horizontal)
sid=info.sv; % Source y position (vertical)
reconSlice=info.vsamp; % Vertical sampling
objDetDistance=info.odd;% Object to detector distance
detLength=info.dpl;% Detector length
numDetPixel=info.dp;% Number of pixels in detector
reconHeight=info.rh;% Actual height of the recon space

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
rows = zeros(2*reconSlice*numProj*numDetPixel,1);
cols = rows;
val1 = rows;
val2 = rows;
idx = 0;

% Initialize the weighting for the measurement
wt = ones (numProj*numDetPixel,1);

% Find the normalization factor
normA = ones(numDetPixel*numProj,2);
for i = 1:numProj
    src = srcX(i);
    for nj = 1:numDetPixel
        dp = (i-1)*numDetPixel+nj;
        botLength = abs(src-(nj-1/2)*pixelSize + detLength/2);
        rBot = (sid - objDetDistance)/sid;
        rTop = (sid - reconHeight - objDetDistance)/sid;
        botPos = botLength*rBot;
        topPos = botLength*rTop;
        
        if abs(src)<= (detLength/2)
            normA(dp,1) = abs(topPos-botPos);
            normA(dp,2) = reconHeight;
        elseif src < (-detLength/2)
            if (botPos <= -(src+detLength/2))
                normA(dp,1) = 0;
                normA(dp,2) = 0;
            elseif (topPos < -(src+detLength/2))
                normA(dp,1) = botPos+(src+detLength/2);
                normA(dp,2) = normA(dp,1)*sid/botLength;
            elseif (topPos >= -(src+detLength/2))
                normA(dp,1) = abs(topPos-botPos);
                normA(dp,2) = reconHeight;
            end
        elseif (src > detLength/2)
            if (botPos <= (src-detLength/2))
                normA(dp,1) = 0;
                normA(dp,2) = 0;
            elseif (topPos < (src-detLength/2))
                normA(dp,1) = botPos-(src-detLength/2);
                normA(dp,2) = normA(dp,1)*sid/botLength;
            elseif (topPos >= (src-detLength/2))
                normA(dp,1) = abs(topPos-botPos);
                normA(dp,2) = reconHeight;
            end
        end
    end
end
normA = (normA/reconSlice).^2;

for i = 1:numProj
    
    % The coordinate of the source location for this projection
    src = srcX(i);

    % Source to detector edge distance
    src_dist = abs(src) - detLength/2;

    % Computing the weighting factor
    if src_dist>0
        % Edge of the recon space rectangular
        top_edge = objDetDistance+reconHeight;
        bot_edge = objDetDistance;

        % The calculated distance should divide by the unit length "pixelSize" to
        % find the pixel location
        det_end = ceil(top_edge/(sid-top_edge)*src_dist/pixelSize);
        det_srt = ceil(bot_edge/(sid-bot_edge)*src_dist/pixelSize);
        ind = det_srt:det_end;
        det_dist = (ind-1/2)*pixelSize;
        edge_height = det_dist./(src_dist + det_dist)*sid;
        w2 = (edge_height - objDetDistance)/reconHeight;

        if src<0
            idx1 = (i-1)*numDetPixel;
            wt(idx1+1:idx1+det_srt-1)=0.001;
            wt(idx1+det_srt:idx1+det_end)=w2;

        elseif src>0
            idx1 = i*numDetPixel;
            wt(idx1-det_srt+2:idx1)=0.001;
            wt(idx1-det_end+1:idx1-det_srt+1)=fliplr(w2);
        end
    end
    wt(wt<0) = 0.001;
    
    %% Start
    if abs(src)<= (detLength/2)
        for j = 1:reconSlice
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
                    val1(idx) = normA(dp,1);
                    val2(idx) = normA(dp,2);
                    cols(idx) = floor(nsp)*reconSlice + j;
                else
                    idx = idx +1;
                    rows(idx) = dp;
                    tmpVal = (ceil(nsp) - nsp)/ratio;
                    val1(idx) = normA(dp,1)*tmpVal^2;
                    val2(idx) = normA(dp,2)*tmpVal^2;
                    cols(idx) = floor(nsp)*reconSlice + j;
                    idx = idx +1;
                    rows(idx) = dp;
                    tmpVal = (nep - floor(nep))/ratio;
                    val1(idx) = normA(dp,1)*tmpVal^2;
                    val2(idx) = normA(dp,2)*tmpVal^2;
                    cols(idx) = floor(nep)*reconSlice + j;
                end
            end
        end
    elseif src<= (-detLength/2)
        for j = 1:reconSlice
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
            dp = (i-1)*numDetPixel+ceil(startPixel);
            rows(idx) = dp;
            tmpVal = startPoint/ratio;
            val1(idx) = normA(dp,1)*tmpVal^2;
            val2(idx) = normA(dp,2)*tmpVal^2;
            cols(idx) = floor(startPixel)*reconSlice + j;
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
                val1(idx) = normA(dp,1);
                val2(idx) = normA(dp,2);
                cols(idx) = floor(nsp)*reconSlice + j;
            else
                idx = idx +1;
                rows(idx) = dp;
                tmpVal = (ceil(nsp) - nsp)/ratio;
                val1(idx) = normA(dp,1)*tmpVal^2;
                val2(idx) = normA(dp,2)*tmpVal^2;
                cols(idx) = floor(nsp)*reconSlice + j;
                idx = idx +1;
                rows(idx) = dp;
                tmpVal = (nep - floor(nep))/ratio;
                val1(idx) = normA(dp,1)*tmpVal^2;
                val2(idx) = normA(dp,2)*tmpVal^2;
                cols(idx) = floor(nep)*reconSlice + j;
            end
        end
        end
    elseif src>= (detLength/2)
        for j = 1:reconSlice
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
                val1(idx) = (endPoint^2)*normA(dp,1);
                val2(idx) = (endPoint^2)*normA(dp,2);
                cols(idx) = floor(endPixel)*reconSlice + j;
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
                    val1(idx) = normA(dp,1);
                    val2(idx) = normA(dp,2);
                    cols(idx) = floor(nsp)*reconSlice + j;
                else
                    idx = idx +1;
                    rows(idx) = dp;
                    tmpVal = (ceil(nsp) - nsp)/ratio;
                    val1(idx) = (tmpVal^2)*normA(dp,1);
                    val2(idx) = (tmpVal^2)*normA(dp,2);
                    cols(idx) = floor(nsp)*reconSlice + j;
                    idx = idx +1;
                    rows(idx) = dp;
                    tmpVal = (nep - floor(nep))/ratio;
                    val1(idx) = (tmpVal^2)*normA(dp,1);
                    val2(idx) = (tmpVal^2)*normA(dp,2);
                    cols(idx) = floor(nep)*reconSlice + j;
                end
            end           
        end
    end
end

% Truncate excess zeros.
rows(idx+1:end) = [];
cols(idx+1:end) = [];
val1(idx+1:end) = [];
val2(idx+1:end) = [];

% Create sparse matrix A from the stored values.
Ax = sparse(rows,cols,val1,numDetPixel*numProj,reconSlice*numDetPixel);
Ay = sparse(rows,cols,val2,numDetPixel*numProj,reconSlice*numDetPixel);