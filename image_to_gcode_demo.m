
imagePath       = 'star12.jpg';    % input image file (put it in current folder)
outputGcode     = 'output.gcode'; % output G-code file

% Physical workspace (physical mm, negative from manual origin)
workWidthPhysMM  = 200;
workHeightPhysMM = 150;

% Desired drawing size in PHYSICAL mm (leave [] to auto-fit)
desiredWidthPhysMM  = [];
desiredHeightPhysMM = [];

% Mechanical scale: G-code mm -> physical mm
% If 1 mm in G-code moves 10 mm physically, use mechanicalScale = 0.1
mechanicalScale = 0.1;

% Image processing + downscale to reduce points
downscaleFactor = 0.6;      % 0 < f <= 1 (set f=1 to keep original)
useAdaptive      = true;    % use imbinarize if available
binThresh        = 0.5;     % fallback threshold
edgeMethod       = 'Canny';
edgeDilatePx     = 0;
minBoundaryPx    = 20;

% Path sampling & simplification
pathStepMM_gcode = 3.0;     % spacing in G-code mm units (bigger -> fewer points)
simplifyPaths    = true;
simplifyTolPx    = 3;       % pixels; fallback decimation if reducepoly not available

% Motion / safety
feedRate        = 200;    % mm/min
rapidSafetyMM   = 2.0;    % small negative park (keeps us away from +0)
insertStops     = true;   % insert M0 before drawing each contour (manual pencil ops)
commentPrefix   = '; ';
fmt             = '%.3f';

% Homing/zeroing flags
assumeManualHome     = true;   % true -> do NOT home, operator must position machine
useG92ToSetZero      = false;  % optional: emit G92 X0 Y0 if you want software zero
returnToHomeAfterRun = false;  % do NOT return to 0,0 at end by default

% Axis inversion (flip signs if your machine moves opposite)
invertX = false;
invertY = false;

% Direction test
writeDirectionTest = true;
safeTestDistance_gcode = 5.0;  % small test move in G-code units (pencil UP)


if ~exist(imagePath,'file')
    error('Input image not found: %s', imagePath);
end

I = imread(imagePath);
if size(I,3) > 1
    I = rgb2gray(I);
end
I = im2double(I);

% optional downscale to reduce points
if exist('downscaleFactor','var') && downscaleFactor>0 && downscaleFactor<1
    I = imresize(I, downscaleFactor, 'bilinear');
end
[imgH, imgW] = size(I);

% Binarize (adaptive if available)
if useAdaptive && exist('imbinarize','file')
    BW1 = imbinarize(I,'adaptive','ForegroundPolarity','dark','Sensitivity',0.45);
    BW2 = imbinarize(I,'adaptive','ForegroundPolarity','bright','Sensitivity',0.45);
    BW1 = bwareaopen(BW1, 8); BW2 = bwareaopen(BW2, 8);
    if nnz(BW2) > nnz(BW1), BW = BW2; else BW = BW1; end
else
    if exist('graythresh','file')
        lvl = graythresh(I);
    else
        lvl = binThresh;
    end
    BW = im2bw(I, lvl); 
    BW = bwareaopen(BW,8);
end
figure, imshow(BW);
% edges
E = edge(BW, edgeMethod);
if edgeDilatePx > 0
    se = strel('disk', edgeDilatePx);
    E = imdilate(E, se);
end
figure, imshow(E);
% boundaries (outer)
[B, ~] = bwboundaries(E, 'noholes');
keep = cellfun(@(b) size(b,1) >= minBoundaryPx, B);
B = B(keep);
if isempty(B)
    error('No contours found. Try simplifying image or adjusting thresholds.');
end


% Fit drawing physically
if isempty(desiredWidthPhysMM) && isempty(desiredHeightPhysMM)
    scaleW = workWidthPhysMM / imgW;
    scaleH = workHeightPhysMM / imgH;
    physScale = min(scaleW, scaleH);
    drawW_phys = imgW * physScale;
    drawH_phys = imgH * physScale;
else
    if isempty(desiredWidthPhysMM)
        drawH_phys = min(desiredHeightPhysMM, workHeightPhysMM);
        drawW_phys = drawH_phys * (imgW / imgH);
        if drawW_phys > workWidthPhysMM, s = workWidthPhysMM/drawW_phys; drawW_phys = drawW_phys*s; drawH_phys = drawH_phys*s; end
    elseif isempty(desiredHeightPhysMM)
        drawW_phys = min(desiredWidthPhysMM, workWidthPhysMM);
        drawH_phys = drawW_phys * (imgH / imgW);
        if drawH_phys > workHeightPhysMM, s = workHeightPhysMM/drawH_phys; drawW_phys = drawW_phys*s; drawH_phys = drawH_phys*s; end
    else
        drawW_phys = min(desiredWidthPhysMM, workWidthPhysMM);
        drawH_phys = min(desiredHeightPhysMM, workHeightPhysMM);
        sW = drawW_phys / imgW; sH = drawH_phys / imgH; physScale = min(sW,sH);
        drawW_phys = imgW * physScale; drawH_phys = imgH * physScale;
    end
end

% convert to G-code units (controller mm)
drawW_gmm = drawW_phys * mechanicalScale;
drawH_gmm = drawH_phys * mechanicalScale;
mmPerPx_gcode = min(drawW_gmm / imgW, drawH_gmm / imgH);

% sanity
if drawW_gmm > workWidthPhysMM * mechanicalScale + 1e-9 || drawH_gmm > workHeightPhysMM * mechanicalScale + 1e-9
    error('Computed G-code extents exceed allowed envelope. Reduce desired size or workspace.');
end

paths = cell(size(B));
for k = 1:numel(B)
    pts_px = B{k}; % Nx2 [row col]
    % ----- inline resamplePolyline (no function def) -----
    % compute segment lengths in pixels
    dP = diff(pts_px,1,1);
    segLen_px = sqrt(sum(dP.^2,2));
    cum_px = [0; cumsum(segLen_px)];
    total_mm = cum_px(end) * mmPerPx_gcode;
    if total_mm < pathStepMM_gcode
        pts_rs = pts_px([1; end], :);
    else
        nSteps = max(2, ceil(total_mm / pathStepMM_gcode) + 1);
        s_mm = linspace(0, total_mm, nSteps);
        s_px = s_mm / mmPerPx_gcode;
        r = interp1(cum_px, pts_px(:,1), s_px, 'linear');
        c = interp1(cum_px, pts_px(:,2), s_px, 'linear');
        pts_rs = [r(:), c(:)];
    end
    % ----- end inline resamplePolyline -----

    % optional simplify: prefer reducepoly if available
    if simplifyPaths && exist('reducepoly','file') == 2
        xy = [pts_rs(:,2), pts_rs(:,1)];  % [x y]
        xy2 = reducepoly(xy, simplifyTolPx);
        pts_rs = [xy2(:,2), xy2(:,1)];
    elseif simplifyPaths
        stride = max(1, round(simplifyTolPx));
        pts_rs = pts_rs(1:stride:end, :);
        if size(pts_rs,1) < 2, pts_rs = pts_px(1:min(2,end),:); end
    end

    r = pts_rs(:,1); c = pts_rs(:,2);
    xg_nom = -((imgW - c) * mmPerPx_gcode); % nominal <= 0
    yg_nom = -(r * mmPerPx_gcode);         % nominal <= 0

    signX = 1; if invertX, signX = -1; end
    signY = 1; if invertY, signY = -1; end
    xg = signX * xg_nom;
    yg = signY * yg_nom;

    if ~invertX, xg = min(xg, 0); end
    if ~invertY, yg = min(yg, 0); end

    paths{k} = [xg(:), yg(:)];
end


% ensure output uses an absolute path in the current folder
outputGcode = fullfile(pwd, outputGcode);  

% try to open and show the system message if it fails
[fid, msg] = fopen(outputGcode, 'w');
if fid == -1
    error('Failed to open "%s" for writing. System message: %s\nCheck your current folder (pwd) and permissions.', outputGcode, msg);
end
fprintf('Writing G-code to: %s\n', outputGcode);


fprintf(fid, '%s Image-to-GCode script (NO Z, NO homing).\n', commentPrefix);
fprintf(fid, '%s Operator must manually set physical corner as X=0 Y=0 before running.\n', commentPrefix);
fprintf(fid, '%s Mechanical scale: 1 mm G-code -> %.6f mm physical\n', commentPrefix, 1/mechanicalScale);
fprintf(fid, 'G21\nG90\n');

if assumeManualHome
    fprintf(fid, '%s ASSUME MANUAL HOME (no homing commands emitted)\n', commentPrefix);
    if useG92ToSetZero
        fprintf(fid, '%s Emitting G92 X0 Y0 (software zero) -- only enable if you want this\n', commentPrefix);
        fprintf(fid, 'G92 X0 Y0\n');
    end
else
    fprintf(fid, '%s WARNING: assumeManualHome=false. Script still avoids homing commands.\n', commentPrefix);
end

% small negative park to avoid rounding to positive values
startX = min(-abs(rapidSafetyMM), 0);
startY = min(-abs(rapidSafetyMM), 0);
if rapidSafetyMM > 0
    fprintf(fid, '%s Move to small negative park to avoid rounding to +0\n', commentPrefix);
    fprintf(fid, 'G0 X%s Y%s\n', num2str(startX, fmt), num2str(startY, fmt));
end

% safe built-in direction test (pencil UP) - does NOT move to +0
if writeDirectionTest
    fprintf(fid, '%s --- DIRECTION TEST (pencil UP) ---\n', commentPrefix);
    fprintf(fid, '%s Operator must be at manual corner (X=0 Y=0) before running this test.\n', commentPrefix);
    fprintf(fid, 'G1 X%s Y%s F%d\n', num2str(-safeTestDistance_gcode, fmt), num2str(startY, fmt), round(feedRate));
    fprintf(fid, 'G1 X%s Y%s F%d\n', num2str(-0.5*safeTestDistance_gcode, fmt), num2str(startY, fmt), round(feedRate));
    fprintf(fid, 'G1 X%s Y%s F%d\n', num2str(-0.5*safeTestDistance_gcode, fmt), num2str(-safeTestDistance_gcode, fmt), round(feedRate));
    fprintf(fid, '%s --- END DIRECTION TEST ---\n', commentPrefix);
end

fprintf(fid, 'F%d\n', round(feedRate));

% write contours
for k = 1:numel(paths)
    P = paths{k};
    if size(P,1) < 2, continue; end

    if insertStops
        fprintf(fid, '%s --- Contour %d : MANUAL pencil lift/lower ---\n', commentPrefix, k);
        fprintf(fid, 'M0\n');
    else
        fprintf(fid, '%s --- Contour %d ---\n', commentPrefix, k);
    end

    fprintf(fid, 'G0 X%s Y%s\n', num2str(P(1,1), fmt), num2str(P(1,2), fmt));

    if insertStops
        fprintf(fid, '%s LOWER pencil now (manual), then resume\n', commentPrefix);
        fprintf(fid, 'M0\n');
    end

    for i = 2:size(P,1)
        fprintf(fid, 'G1 X%s Y%s\n', num2str(P(i,1), fmt), num2str(P(i,2), fmt));
    end

    if rapidSafetyMM > 0
        fprintf(fid, 'G0 X%s Y%s\n', num2str(min(P(end,1), startX), fmt), num2str(min(P(end,2), startY), fmt));
    end
end

% final park (do not return to 0,0)
if rapidSafetyMM > 0
    fprintf(fid, '%s Final small negative park\n', commentPrefix);
    fprintf(fid, 'G0 X%s Y%s\n', num2str(startX, fmt), num2str(startY, fmt));
end
if returnToHomeAfterRun
    fprintf(fid, '%s Returning to X0 Y0 at end (ONLY set true if safe)\n', commentPrefix);
    fprintf(fid, 'G0 X0 Y0\n');
else
    fprintf(fid, '%s NOT returning to X0 Y0 (operator must park manually)\n', commentPrefix);
end
fprintf(fid, 'M2\n');
fclose(fid);
% show small preview so you can confirm it was written
fprintf('Wrote G-code to: %s\n', outputGcode);
try
    fprintf('--- G-code preview (first 40 lines) ---\n');
    fid2 = fopen(outputGcode,'r');
    for i=1:40
        t = fgetl(fid2);
        if ~ischar(t), break; end
        disp(t);
    end
    fclose(fid2);
catch
    % ignore preview errors
end


fprintf('G-code written to: %s\n', outputGcode);
fprintf('Drawing physical size: W=%.1f mm, H=%.1f mm\n', drawW_phys, drawH_phys);
fprintf('Drawing G-code units size: W=%.1f mm, H=%.1f mm\n', drawW_gmm, drawH_gmm);


