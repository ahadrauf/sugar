function writebeam(fid, Rp, R, L, W, H, var_ids)

fwrite(fid, 1, 'integer*4');  % Identify display model type
fwrite(fid, R', 'real*4');    % Output rotation transposed (column-major)
fwrite(fid, L, 'real*4');     % Length
fwrite(fid, W, 'real*4');     % Width
fwrite(fid, H, 'real*4');     % Height
fwrite(fid, Rp, 'real*4');    % "Root" point
fwrite(fid, var_ids-1, 'integer*4');  % Variable indices

