function dem = viz_tif(tifFile, varargin)
%VIZ_TIF Load a DEM GeoTIFF and optionally visualize it.
%   dem = viz_tif(tifFile) loads tifFile, removes NoData, and returns a
%   struct with downsampled X/Y/Z grids, bounds, and a griddedInterpolant.
%
%   dem = viz_tif(tifFile, 'sampleStep', N, 'pointCloudStep', M,
%   'makePlots', true) adjusts the raster decimation used for both the
%   returned grid and optional plots. Use sampleStep=1 for full fidelity.

p = inputParser;
addParameter(p, 'sampleStep', 10, @(v) isnumeric(v) && isscalar(v) && v >= 1);
addParameter(p, 'pointCloudStep', 500, @(v) isnumeric(v) && isscalar(v) && v >= 1);
addParameter(p, 'makePlots', false, @(v) islogical(v) && isscalar(v));
parse(p, varargin{:});
opts = p.Results;

[Zfull, R] = readgeoraster(tifFile);
Zfull = double(Zfull);

% Handle NoData (if present)
info = geotiffinfo(tifFile);
nodataVal = [];
if isfield(info, 'GeoTIFFTags') && isfield(info.GeoTIFFTags, 'GDAL_NODATA')
    nodataVal = str2double(info.GeoTIFFTags.GDAL_NODATA);
    Zfull(Zfull == nodataVal) = NaN;
end

Zfull(Zfull < 0) = NaN;   % just in case
Zfull = fillmissing(Zfull, 'nearest'); % avoid holes that produce NaNs later

% Build X/Y coordinate grid (UTM meters)
x = R.XLimWorld(1) + (0:R.RasterSize(2)-1) * R.CellExtentInWorldX;
y = R.YLimWorld(2) - (0:R.RasterSize(1)-1) * R.CellExtentInWorldY;
[Xfull, Yfull] = meshgrid(x, y);

% Ensure ascending coordinates for griddedInterpolant (flip if descending)
if numel(x) > 1 && x(2) < x(1)
    x = fliplr(x);
    Xfull = fliplr(Xfull);
    Zfull = fliplr(Zfull);
end
if numel(y) > 1 && y(2) < y(1)
    y = flipud(y(:))';
    Yfull = flipud(Yfull);
    Zfull = flipud(Zfull);
end

% Downsample for speed / accuracy studies
step = max(1, round(opts.sampleStep));
X = Xfull(1:step:end, 1:step:end);
Y = Yfull(1:step:end, 1:step:end);
Z = Zfull(1:step:end, 1:step:end);

dem.X = X;
dem.Y = Y;
dem.Z = Z;
dem.R = R;
dem.bounds.xlim = [min(X(:)), max(X(:))];
dem.bounds.ylim = [min(Y(:)), max(Y(:))];
dem.bounds.center = [mean(dem.bounds.xlim), mean(dem.bounds.ylim)];
dem.sampleStep = step;
dem.nodataVal = nodataVal;
dem.interp = griddedInterpolant(X', Y', Z', "linear", "nearest");

if opts.makePlots
    figure;
    surf(X, Y, Z, Z, 'EdgeColor', 'none');  % color = elevation
    colormap(parula);
    cb = colorbar;
    cb.Label.String = 'Elevation (m)';
    cb.Label.FontSize = 12;
    title('3D DEM Surface');
    xlabel('Easting (m, UTM Zone 12N)');
    ylabel('Northing (m, UTM Zone 12N)');
    zlabel('Elevation (m)');
    axis equal;
    view(45, 60);
    camlight headlight;

    figure;
    valid = ~isnan(Zfull);
    pcX = Xfull(valid);
    pcY = Yfull(valid);
    pcZ = Zfull(valid);
    pcStep = max(1, round(opts.pointCloudStep));
    pcX = pcX(1:pcStep:end);
    pcY = pcY(1:pcStep:end);
    pcZ = pcZ(1:pcStep:end);
    scatter3(pcX, pcY, pcZ, 5, pcZ, 'filled');  % colored by elevation
    colormap(parula);
    cb2 = colorbar;
    cb2.Label.String = 'Elevation (m)';
    xlabel('Easting (m, UTM Zone 12N)');
    ylabel('Northing (m, UTM Zone 12N)');
    zlabel('Elevation (m)');
    title('DEM Point Cloud');
    axis equal;
    view(3);
    grid on;
end

end
