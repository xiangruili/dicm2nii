function varargout = nii_viewer(fname, varargin)
% Basic tool to visualize NIfTI images.
% 
%  NII_VIEWER('/data/subj2/fileName.nii.gz')
%  NII_VIEWER('background.nii', 'overlay.nii')
%  NII_VIEWER('background.nii', {'overlay1.nii' 'overlay2.nii'})
% 
% If no input is provided, the viewer will load included MNI_2mm brain as
% background NIfTI. Although the preferred format is NIfTI, NII_VIEWER accepts
% any files that can be converted into NIfTI by dicm2nii, including NIfTI,
% dicom, PAR/REC, etc. In case of CIfTI file, it will show both the volume and 
% surface view if GIfTI is available.
% 
% Here are some features and usage.
% 
% The basic use is to open a NIfTI file to view. When a NIfTI (background) is
% open, the display always uses the image plane close to xyz axes (voxel space)
% even for oblique acquisition. The possible confusion arises if the acquisition
% was tilted with a large angle, and then the orientation labels will be less
% accurate. The benefit is that no interpolation is needed for the background
% image, and there is no need to switch coordinate system when images with
% different systems are overlaid. The display is always in correct scale at
% three vies even with non-isotropic voxels. The displayed IJK always correspond
% to left -> right, posterior -> anterior and inferior -> superior directions,
% while the NIfTI data may not be saved in this order or along these directions.
% The I-index is increasing from left to right even when the display is flipped
% as radiological convention (right on left side).
% 
% Navigation in 4D can be by mouse click, dialing IJK and volume numbers, or
% using keys (arrow keys and [ ] for 3D, and < > for volume).
% 
% After the viewer is open, one can drag-and-drop a NIfTI file into the viewer
% to open as background, or drop a NIfTI with Control key down to add it as
% overlay.
% 
% By default, the viewer shows full view of the background image data. The
% zoom-in always applies to three views together, and enlarges around the
% location of the crosshair. To zoom around a different location, set the
% crosshair to the interested location, and apply zoom again either by View ->
% Zoom in, or pressing Ctrl (Cmd) and +/-. View -> Set crosshair at -> center of
% view (or pressing key C) can set the crosshair to the center of display for
% all three views.
% 
% Overlays are always mapped onto the coordinate of background image, so
% interpolation (nearest/linear/cubic/spline) is usually involved. The overlay
% makes sense only when it has the same coordinate system as the background
% image, while the resolution and dimension can be different. The viewer tries
% to match any of sform and qform between the images. If there is no match, a
% warning message will show up.
% 
% A special overlay feature "Add aligned overlay" can be used to check the
% effect of registration, or overlay an image to a different coordinate system
% without creating a transformed file. It will ask for two files. The first is
% the overlay NIfTI file, and the second is either a transformation matrix file
% which aligns the overlay to the background image, or a warp file which
% transforms the overlay into the background reference.
% 
% Here is an example to check FSL alignment. From a .feat/reg folder, Open
% "highres" as background image. "Add overlay" for "example_func". If there is
% head movement between the highres and the functional scans, the overlap will
% be off. Now "Add aligned overlay" for "example_func", and use
% example_func2highres.mat as the matrix. The two dataset should overlap well if
% the alignment matrix is accurate.
% 
% Here is another example to overlay to a different coordinate system for FSL
% result. Open .feat/reg/standard.nii.gz as background. If an image like
% .feat/stats/zstat1.nii.gz is added as an overlay, a warning message will say
% inconsistent coordinate system since the zstat image is in Scanner Anat. The
% correct way is to "Add aligned overlay" for zstat1, and either use
% .feat/reg/example_func2standard.mat for linear transformation or better use
% .feat/reg/example_func2standard_warp.nii.gz if available for alignment.
% 
% When the mouse pointer is moved onto a voxel, the x/y/z coordinates and voxel
% value will show on the panel. If there is an overlay, the overlay voxel value
% will also show up, unless its display is off. When the pointer moves onto the
% panel or out of an image, the information for the voxel at crosshair will be
% displayed. The display format is as following with val of the top image
% displayed first:
%  (x,y,z): val_1 val_2 val_3 ...
%
% Note that although the x/y/z coordinates are shared by background and overlay
% images, IJK indices are always for background image (name shown on title bar).
% 
% The mouse-over display can be turned on/off from Help -> Preferences ...
% 
% If there is a .txt label file in the same folder as the NIfTI file, like for
% AAL, the labels will be shown instead of voxel value. The txt file should have
% a voxel value and ROI label pair per line, like
%  1 Precentral_L
%  2 Precentral_R
%  3 ...
% 
% Image display can be smoothed in 3D (default is off). The smooth is slow when
% the image dimension is large, even when the current implementation of smooth
% does not consider voxel size.
% 
% Background image and overlays are listed at the top-left of the panel. All
% parameters of the bottom row of the panel are for the highlighted image. This
% feature is indicated by the blue background of the parameters. Each NIfTI file
% has its own set of parameters (display min and max value, LUT, alpha, whether
% to smooth, interpolation method, and volume number) to control its display.
% Moving the mouse onto a parameter will show the meaning of the parameter.
% 
% The lastly added overlay is on the top of display and top of the file list.
% This also applies in case there is surface figure for CIfTI files. The
% background and overlay order can be changed by the two small arrows next to
% the list, or from Overlay -> Move selected image ...
% 
% Each NIfTI display can be turned on/off by clicking the small checkbox at the
% left side of the file (or pressing spacebar for the selected NIfTI). This
% provides a way to turn on/off an overlay to view the overlap. Most operations
% are applied to the selected NIfTI in the list, such as Show NIfTI hdr/ext
% under Window menu, Move/Close overlay under Overlay menu, and most operations
% under File menu.
% 
% A NIfTI mask can be applied to the selected image. Ideally, the mask should be
% binary, and the image corresponding to the non-zero part of the mask will be
% displayed. If non-binary mask is detected, a threshold to binarize will be
% asked. If the effect is not satisfied with the threshold, one can apply the
% mask with a different threshold without re-loading image. The way to remove a
% mask is to Remove, then Add overlay again. In case one likes to modulate the
% selected image with another NIfTI image (multiply two images), File -> Apply
% modulation will do it. If the mask image is not within 0 and 1, the lower and
% upper bound will be asked to normalize the range to [0 1]. A practical use of
% modulation is to use dti_FA map as weight for dti_V1 RGB image.
% 
% For multi-volume data, one can change the Volume Number (the parameter at
% rightmost of the panel) to check the head motion. Click in the number dialer
% or press < or > key, to simulate movie play. It is better to open the 4D data
% as background, since it may be slower to map it to the background image.
% 
% Popular LUT options are implemented. Custom LUT can be added by Overlay ->
% Load LUT for selected overlay. The custom LUT file can be in text format (each
% line represents an RGB triplet, while the first line corresponds to the value
% 0 in the image data), or binary format (uint8 for all red, then green then
% blue). The color coding can be shown by View -> Show colorbar. There are
% several special LUTs. The "two-sided" allows to show both positive and
% negative data in one view. For example, if the display range is 3 to 10 for a
% t-map, positive T above 3 will be coded as red-yellow, and T below -3 will be
% coded as blue-green. This means the absolute display range values are used.
% 
% One of the special LUT is "lines". This is for normalized vector display,
% e.g. for diffusion vector. The viewer will refuse the LUT selection if
% the data is not normalized vector. Under this LUT, all other parameters
% for the display are ignored. The color of the "lines" is the max color of
% previous LUT. For example, if one likes to show blue vector lines, choose
% LUT "blue" first, then change it to "lines".
% 
% In case of complex image, most LUTs will display only the magnitude of the
% data, except the following three phase LUTs, where the magnitude is used as
% mask. Here is how the 3 phase LUTs encodes phase from 0 to 360 degree:
%  phase:  red-yellow monotonically,
%  phase3: red-yellow-green-yellow-red circularly, and 
%  phase6: red-yellow-green/violet-blue-cyan, with sharp change at 180 degree. 
% These LUTs are useful to visualize activation of periodic stimulus, such as
% those by expanding/shrinking ring or rotating wedge for retinotopy. To use
% this feature, one can save an complex NIfTI storing the Fourier component at
% the stimulus frequency.
% 
% Note that, for RGB NIfTI, the viewer always displays the data as color
% regardless of the LUT option. The display min and max value also have no
% effect on RGB image. There is a special LUT, RGB, which is designed to display
% non-RGB NIfTI data as RGB, e.g. FSL-style 3-volome data. 
% 
% The viewer figure can be copied into clipboard (not available for Linux) or
% saved as variety of image format. For high quality picture, one can increase
% the output resolution by Help -> Preferences -> Resolution. Higher resolution
% will take longer time to copy or save, and result in larger file. If needed,
% one can change to white background for picture output. With white background,
% the threshold for the background image needs to be carefully selected to avoid
% strange effect with some background images. For this reason, white background
% is intended only for picture output.
% 
% The selected NIfTI can also be saved as different format from File -> Save
% NIfTI as. For example, a file can be saved as a different resolution. With a
% transformation matrix, a file can also be saved into a different template. The
% latter is needed for FSLview since it won't allow overlay with different
% resolution or dimension at least till version 5.0.8.
% 
% See also NII_TOOL DICM2NII NII_XFORM
 
%% By Xiangrui Li (xiangrui.li at gmail.com)
% 151021 Almost ready to include into dicm2nii package.
% 191019 start to use uifigure.
% 260129 uifigure version replaces figure version.
%%
 
if nargin==2 && ischar(fname) && strcmp(fname, 'func_handle')
    varargout{1} = str2func(varargin{1});
    return;
elseif nargin>1 && ischar(fname) && strcmp(fname, 'LocalFunc')
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
 
if nargin<1 || isempty(fname) % open the included standard_2mm
    fname = fullfile(fileparts(mfilename('fullpath')), 'example_data.mat'); 
    fname = load(fname, 'nii'); fname = fname.nii;
end
 
nii = get_nii(fname);
[p, hs.form_code, rg, dim] = read_nii(nii); % re-oriented
p.Ri = inv(p.R);
nVol = size(p.nii.img, 4);
hs.bg.R   = p.R;
hs.bg.Ri  = p.Ri;
hs.bg.hdr = p.hdr0;
if ~isreal(p.nii.img)
    p.phase = angle(p.nii.img); % -pi to pi
    p.phase = mod(p.phase/(2*pi), 1); % 0~1
    p.nii.img = abs(p.nii.img); % real now
end
hs.dim = single(dim); % single saves memory for ndgrid
hs.pixdim = p.pixdim;
hs.gap = min(hs.pixdim) ./ hs.pixdim * 3; % in unit of smallest pixdim
 
p.lb = rg(1); p.ub = rg(2);
p = dispPara(p);
[pName, niiName, ext] = fileparts(p.nii.hdr.file_name);
if strcmpi(ext, '.gz'), [~, niiName] = fileparts(niiName); end
 
if nargin>1 && any(ishandle(varargin{1})) % called by Open or dnd
    fh = varargin{1};
    hsN = guidata(fh);
    pf = hsN.pref.UserData; % use the current pref unless for new figure
    close(fh);
else
    pf = getpref('nii_viewer_para');
    if isempty(pf) || ~isfield(pf, 'layout') % check lastly-added field
        pf = struct('openPath', pwd, 'addPath', pwd, 'interp', 'linear', ...
            'extraV', NaN, 'dpi', '0', 'rightOnLeft', false, ...
            'mouseOver', true, 'layout', 1);
        setpref('nii_viewer_para', fieldnames(pf), struct2cell(pf));
    end
end
[siz, axPos, figPos] = plot_pos(dim.*hs.pixdim, pf.layout);

figNam = char(p.nii.hdr.file_name);
if numel(figNam)>40, figNam = [figNam(1:20) '...' figNam((-19:0)+end)]; end
figNam = ['nii_viewer - ' figNam ' (' formcode2str(hs.form_code(1)) ')'];
fh = uifigure('Name', figNam, 'Position', [figPos siz+[2 70]]);
c = round(p.Ri * [0 0 0 1]'); c = c(1:3)' + 1; % 
ind = c<=1 | c>=dim;
c(ind) = round(dim(ind)/2);
% c = round(dim/2); % start cursor at the center of images
xyz = round(p.R * [c-1 1]'); % take care of rounding error
hs.fig = fh;
cb = @(cmd) {@nii_viewer_cb cmd fh}; % callback shortcut
if nargout, varargout{1} = fh; end
 
%% split into control panel and image frame
L0 = uigridlayout(fh, [2 1], 'RowHeight', {68 '1x'}, 'Padding', 1, 'RowSpacing', 0);
hs.panel = uipanel(L0, 'BorderType', 'none'); % control panel
hs.frame = uipanel(L0, 'BorderType', 'none', 'BackgroundColor', 'k', 'UserData', siz);
try uiFileDnD(hs.frame, cb("drop")); catch me, disp(me.message); end
try fh.Icon = fullfile(fileparts(mfilename('fullpath')), 'icon.png'); end % since R2020b
 
% Nx2 uitable simulate checkboxList
L12 = uigridlayout(hs.panel, [1 2], 'ColumnWidth', {108 '1x'}, 'Padding', 0, 'ColumnSpacing', 0);
hs.files = uitable(L12, 'Data', {true niiName}, 'ColumnEditable', [true false], ...
    'FontName', 'Tahoma', 'FontSize', 11, 'RowName', {}, ...
    'ColumnName', {}, 'ColumnWidth', {20 84}, 'RowStriping', 'off', ...
    'Tooltip', {'Select image to show/modify its display parameters.' ...
        'Click checkbox to turn on/off image'}, ...
    'CellSelectionCallback', cb("file"), 'CellEditCallback', cb("toggle"));
try hs.files.Selection = 1; end
stl = uistyle('BackgroundColor', [124 169 195]/256, 'FontColor', 'w');
addStyle(hs.files, stl, 'row', 1);
 
% panel for controls except uitable
L21 = uigridlayout(L12, [2 1], 'Padding', 0, 'RowSpacing', 2);
L1 = uigridlayout(L21, [1 8], 'Padding', [0 6 6 6], 'ColumnSpacing', 2, ...
    'ColumnWidth', {24 14 50 14 50 14 50 '1x'});
try % since R2020b
    L2 = uigridlayout(L21, [1 8], 'BackgroundColor', stl.BackgroundColor);
catch
    p0 = uipanel(L21, 'BackgroundColor', stl.BackgroundColor); % color purpose only
    L2 = uigridlayout(p0, [1 8]);
end
set(L2, 'Padding', 6, 'ColumnSpacing', 2, 'ColumnWidth', {60 60 90 48 60 76 54 '1x'});
 
hs.stack = uispinner(L1, 'Value', -1, 'HorizontalAlignment', 'left', ...
    'ValueChangedFcn', cb("stack"), 'Limits', [-2 -1], ...
    'Enable', false, 'Tooltip', 'Move highlighted image up/down');
 
labls = 'IJK';
str = ["Left to Right" "Posterior to Anterior" "Inferior to Superior"];
for i = 1:3
    txt = str{i}+", 1:"+dim(i);
    uilabel(L1, 'Text', labls(i), 'HorizontalAlignment', 'right', ...
        'Tooltip', txt, 'FontWeight', 'bold');
    hs.ijk(i) = uispinner(L1, 'Value', c(i), 'ValueDisplayFormat', '%i', ...
        'Limits', [1 max(dim(i),1+1e-6)], 'ValueChangedFcn', cb("ijk"), ...
        'Tooltip', txt, 'FontSize', 11, 'RoundFractionalValues', 'on');
end
hs.value = uilabel(L1, 'Text', '', 'FontSize', 11, ...
    'HorizontalAlignment', 'center', 'Tooltip', '(x,y,z): top ... bottom');
 
% Controls for each file
hs.lb = uispinner(L2, 'Value', p.lb, 'Step', p.lb_step, 'ValueDisplayFormat', '%.4g', ...
    'ValueChangedFcn', cb("lb"), 'Tooltip', 'min value (threshold)');
hs.ub = uispinner(L2, 'Value', p.ub, 'Step', p.ub_step, 'ValueDisplayFormat', '%.5g', ...
    'ValueChangedFcn', cb("ub"), 'Tooltip', 'max value (clipped)');
 
hs.lut = uidropdown(L2, 'Items', p.lut, 'Value', p.lut, ...
    'ValueChangedFcn', cb("lut"), 'Tooltip', 'Lookup table options');
hs.lut.UserData.lastLut = "red";
hs.lut.UserData.lutStr = ["grayscale" "red" "green" "blue" "violet" "yellow" "cyan" ...
    "red-yellow" "blue-green" "two-sided" ...
    "parula" "jet" "hsv" "hot" "cool" "spring" "summer" "autumn" "winter" ...
    "bone" "copper" "pink" "prism" "flag" "phase" "phase3" "phase6" "RGB" "lines" "custom"];
if p.lut=="custom" || size(p.nii.img,8)>1, hs.lut.Enable = 'off'; end
 
hs.alpha = uispinner(L2, 'Value', 1, 'Step', 0.1, 'Limits', [0 1], 'ValueDisplayFormat', '%.1g', ...
    'ValueChangedFcn', cb("alpha"), 'Tooltip', 'Alpha: 0 transparent, 1 opaque');
 
hs.smooth = uicheckbox(L2, 'Value', p.smooth, 'FontColor', 'w', 'Text', 'smooth', ...
    'ValueChangedFcn', cb("smooth"), 'Tooltip', 'Smooth image in 3D');
hs.interp = uidropdown(L2, 'Value', p.interp, 'Enable', 'off', ...
    'Items', ["nearest" "linear" "cubic" "spline"], 'ValueChangedFcn', cb("interp"), ... 
    'Tooltip', 'Interpolation algorithm for overlay');
hs.volume = uispinner(L2, 'Limits', [1 max(nVol,1+1e-6)], 'ValueDisplayFormat', '%i', ...
    'ValueChangedFcn', cb("volume"), 'RoundFractionalValues', 'on', ...
    'Enable', nVol>1, 'Tooltip', ['Volume number, 1:' num2str(nVol)]);
 
%% Three views: sag, cor, tra
drawnow; % update Position for hs.frame/hs.panel
for i = 1:3
    j = 1:3; j(j==i) = [];
    hs.ax(i) = axes(hs.frame, 'Position', axPos(i,:));
    hs.hsI(i) = image(hs.ax(i), zeros(dim(j([2 1])), 'single'));
    hs.ax(i).Toolbar.Visible = 'off';
    hs.ax(i).DataAspectRatio = [1./hs.pixdim(j) 1];
    hold(hs.ax(i), 'on');
    
    x = [c(j(1))+[-1 1 0 0]*hs.gap(j(1)); 0 dim(j(1))+1 c(j(1))*[1 1]];
    y = [c(j(2))+[0 0 -1 1]*hs.gap(j(2)); c(j(2))*[1 1] 0 dim(j(2))+1];
    hs.cross(i,:) = line(hs.ax(i), x, y);
 
    hs.xyz(i) = text(hs.ax(i),0.02, 0.95, num2str(xyz(i)), ...
        'Units', 'normalized', 'FontSize', 12);
end
set(hs.hsI, 'ButtonDownFcn', cb("mousedown"));
p.hsI = hs.hsI; % background img
p.hsI(1).UserData = p; % store everything in sag img UserData
 
labls = 'ASLSLP'; 
pos = [0.96 0.5; 0.47 0.95; 0 0.5; 0.47 0.95; 0 0.5; 0.47 0.05]; 
for i = 1:numel(labls)
    hs.ras(i) = text(hs.ax(ceil(i/2)), pos(i,1), pos(i,2), labls(i), ...
        'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
end
 
hs.ax(4) = axes(hs.frame, 'Position', axPos(4,:));
hs.colorbar = colorbar(hs.ax(4), 'Color', [1 1 1], ...
    'Location', 'West', 'PickableParts', 'none', 'Visible', 'off');
 
% image() reverses YDir. Turn off ax and ticks
set(hs.ax, 'YDir', 'normal', 'Visible', 'off');
set([hs.ras hs.cross(:)' hs.xyz], 'Color', 'b'); %, 'UIContextMenu', '');
set([hs.ras hs.cross(:)' hs.xyz], 'PickableParts', 'none');
 
%% menus
h = uimenu(fh, 'Label', '&File');
uimenu(h, 'Label', 'Open', 'Accelerator', 'O', 'UserData', pName, 'Callback', cb("open"));
uimenu(h, 'Label', 'Open in new window', 'Callback', cb("open"));
uimenu(h, 'Label', 'Apply mask', 'Callback', @addMask);
uimenu(h, 'Label', 'Apply modulation', 'Callback', @addMask);
h_savefig = uimenu(h, 'Label', 'Save figure as');
h_saveas = uimenu(h, 'Label', 'Save NIfTI as');
uimenu(h, 'Label', 'Save volume as ...', 'Callback', cb("saveVolume"));
uimenu(h, 'Label', 'Export as movie ...', 'Callback', cb('MP4'));
uimenu(h, 'Label', 'Crop below crosshair', 'Callback', cb("cropNeck"));
uimenu(h, 'Label', 'Create ROI file ...', 'Callback', cb("ROI"));
uimenu(h, 'Label', 'Close window', 'Accelerator', 'W', 'Callback', 'close gcf');
 
uimenu(h_saveas, 'Label', 'SPM 3D NIfTI (one file/pair per volume)', 'Callback', @save_nii_as);
uimenu(h_saveas, 'Label', 'NIfTI standard RGB (for AFNI, later mricron)', ...
    'Callback', @save_nii_as, 'Separator', 'on');
uimenu(h_saveas, 'Label', 'FSL style RGB (RGB saved in dim 4)', 'Callback', @save_nii_as);
uimenu(h_saveas, 'Label', 'Old mricron style RGB (RGB saved in dim 3)', 'Callback', @save_nii_as);
uimenu(h_saveas, 'Label', 'a copy', 'Callback', @save_nii_as, 'Separator', 'on');
uimenu(h_saveas, 'Label', 'file with a new resolution', 'Callback', @save_nii_as);
uimenu(h_saveas, 'Label', 'file matching background', 'Callback', @save_nii_as);
uimenu(h_saveas, 'Label', 'file in aligned template space', 'Callback', @save_nii_as);
 
fmt = {'png' 'jpg' 'tif' 'bmp' 'pdf' 'eps'};
if ispc, fmt = [fmt 'emf']; end
for i = 1:numel(fmt)
    uimenu(h_savefig, 'Label', fmt{i}, 'Callback', cb("save"));
end
 
h = uimenu(fh, 'Label', '&Edit');
uimenu(h, 'Label', 'Copy figure', 'Callback', cb("copy"));
 
h_over = uimenu(fh, 'Label', '&Overlay');
uimenu(h_over, 'Label', 'Add overlay', 'Accelerator', 'A', 'Callback', cb("add"));
uimenu(h_over, 'Label', 'Add aligned overlay', 'Callback', cb("add"));
hs.overlay(4) = uimenu(h_over, 'Label', 'Manual alignment', 'Callback', ...
    @manual_align, 'Enable', 'off');
 
h = uimenu(h_over, 'Label', 'Move selected image', 'Enable', 'off');
uimenu(h, 'Label', 'to top',         'Callback', cb("stack"), 'Tag', 'top');
uimenu(h, 'Label', 'to bottom',      'Callback', cb("stack"), 'Tag', 'bottom');
uimenu(h, 'Label', 'one level up',   'Callback', cb("stack"), 'Tag', 'up');
uimenu(h, 'Label', 'one level down', 'Callback', cb("stack"), 'Tag', 'down');
hs.overlay(1) = h;
 
hs.overlay(5) = uimenu(h_over, 'Label', 'Remove overlay', 'Accelerator', 'R', ...
    'Callback', cb("close"), 'Enable', 'off');
hs.overlay(2) = uimenu(h_over, 'Label', 'Remove overlays', ...
    'Callback', cb("closeAll"), 'Enable', 'off');
hs.overlay(3) = uimenu(h_over, 'Label', 'Load LUT for current overlay', 'Callback', cb("custom"));
 
h_view = uimenu(fh, 'Label', '&View');
h = uimenu(h_view, 'Label', 'Zoom in by');
for i = [1 1.2 1.5 2 3 4 5 8 10 20]
    uimenu(h, 'Label', num2str(i), 'Callback', cb("zoom"));
end
h = uimenu(h_view, 'Label', 'Layout', 'UserData', pf.layout);
uimenu(h, 'Label', 'one-row', 'Callback', cb("layout"), 'Tag', '1');
uimenu(h, 'Label', 'two-row sag on right', 'Callback', cb("layout"), 'Tag', '2');
uimenu(h, 'Label', 'two-row sag on left', 'Callback', cb("layout"), 'Tag', '3');
uimenu(h_view, 'Label', 'White background', 'Callback', cb("background"));
hLR = uimenu(h_view, 'Label', 'Right on left side', 'Callback', cb("flipLR"));
uimenu(h_view, 'Label', 'Show colorbar', 'Callback', cb("colorbar"));
uimenu(h_view, 'Label', 'Show crosshair', 'Separator', 'on', ...
    'Checked', 'on', 'Callback', cb("cross"));
h = uimenu(h_view, 'Label', 'Set crosshair at');
uimenu(h, 'Label', 'center of view', 'Callback', cb("viewCenter"));
uimenu(h, 'Label', 'center of image', 'Callback', cb("center"));
uimenu(h, 'Label', 'COG of image', 'Callback', cb("cog"));
uimenu(h, 'Label', 'Smoothed maximum', 'Callback', cb("maximum"));
uimenu(h, 'Label', 'a point [x y z] ...', 'Callback', cb("toXYZ"));
uimenu(h, 'Label', 'a point with value of ...', 'Callback', cb("toValue"));
uimenu(h_view, 'Label', 'Crosshair color', 'Callback', cb("color"));
h = uimenu(h_view, 'Label', 'Crosshair gap');
for i = [0 1 2 3 4 5 6 8 10 20 40]
    str = num2str(i); if i==6, str = [str ' (default)']; end %#ok
    uimenu(h, 'Label', str, 'Callback', cb("gap"));
end
h = uimenu(h_view, 'Label', 'Crosshair thickness');
uimenu(h, 'Label', '0.5 (default)', 'Callback', cb("thickness"));
for i = [0.75 1 2 4 8]
    uimenu(h, 'Label', num2str(i), 'Callback', cb("thickness"));
end
 
h = uimenu(fh, 'Label', '&Window');
uimenu(h, 'Label', 'Show NIfTI essentials', 'Callback', cb("essential"));
uimenu(h, 'Label', 'Show NIfTI hdr', 'Callback', cb("hdr"));
uimenu(h, 'Label', 'Show NIfTI ext', 'Callback', cb("ext"));
uimenu(h, 'Label', 'DICOM to NIfTI converter', 'Callback', 'dicm2nii', 'Separator', 'on');
th = uimenu(h, 'Label', 'Time course ...', 'Callback', cb("tc"), 'Separator', 'on');
setappdata(th, 'radius', 6);
th = uimenu(h, 'Label', 'Standard deviation ...', 'Callback', cb("tc"));
setappdata(th, 'radius', 6);
uimenu(h, 'Label', 'Histogram', 'Callback', cb("hist"));
 
h = uimenu(fh, 'Label', '&Help');
hs.pref = uimenu(h, 'Label', 'Preferences', 'UserData', pf, 'Callback', @pref_dialog);
uimenu(h, 'Label', 'Key shortcut', 'Callback', cb("keyHelp"));
uimenu(h, 'Label', 'Show help text', 'Callback', 'doc nii_viewer');
checkUpdate = dicm2nii('', 'checkUpdate', 'func_handle');
uimenu(h, 'Label', 'Check update', 'Callback', @(~,~)checkUpdate('nii_viewer'));
uimenu(h, 'Label', 'About', 'Callback', cb("about"));

%% finalize gui
guidata(fh, hs); % store handles and data
table_width(hs.files);
set(fh, 'WindowKeyPressFcn', @KeyPressFcn, 'CloseRequestFcn', cb("closeFig"));
 
nii_viewer_cb(hs.files, struct('Indices', [1 2]), "file", fh);
if pf.mouseOver, fh.WindowButtonMotionFcn = cb("mousemove"); end
if pf.rightOnLeft, nii_viewer_cb(hLR, [], "flipLR", fh); end
set_cdata(hs);
set_xyz(hs);
 
if nargin>1
    in1 = varargin{1};
    if ischar(in1) || isStringScalar(in1) || isstruct(in1)
        addOverlay(in1, fh);
    elseif iscell(in1) || isstring(in1)
        for i=1:numel(in1), addOverlay(in1{i}, fh); end
    end
end
 
if hs.form_code(1)<1
    uialert(fh, ['There is no valid form code in NIfTI. The orientation ' ...
        'labeling is likely meaningless.'], 'Invalid Form Code');
end
if isfield(p.nii, 'cii'), cii_view(hs); end
 
%% Get info from sag img UserData
function p = get_para(hs, iFile)
if nargin<2, iFile = -hs.stack.Value; end
hsI = findobj(hs.ax(1), 'Type', 'image', '-or', 'Type', 'quiver');
p = hsI(iFile).UserData;
 
%% callbacks
function nii_viewer_cb(h, evt, cmd, fh)
hs = guidata(fh);
if cmd == "ijk" % IJK spinner
    if isempty(h), iAx = 1:3; else, iAx = find(h == hs.ijk); end
    for ix = iAx
        set_cdata(hs, ix);
        set_cross(hs, ix);
    end
    xyz = set_xyz(hs);
    for i = 1:3, hs.xyz(i).String = xyz(i); end % need 3 if oblique
elseif cmd == "mousedown" % image clicked
    ax = hs.fig.CurrentAxes;
    c = round(ax.CurrentPoint(1, 1:2));
    i = 1:3;
    i(ax==hs.ax(1:3)) = [];
    try hs.ijk(i(1)).Value = c(1); catch, return; end
    nii_viewer_cb(hs.ijk(i(1)), [], "ijk", fh)
    hs.ijk(i(2)).Value = c(2);
    nii_viewer_cb(hs.ijk(i(2)), [], "ijk", fh)
elseif any(cmd == ["lb" "ub" "lut" "alpha" "smooth" "interp" "volume"])
    p = get_para(hs);
    val = h.Value;
    if ischar(val), val = string(val); end
    
    if cmd == "smooth" && val && numel(p.nii.img(:,:,:,1))<2
        h.Value = false; return;
    elseif cmd == "lut" && val == "lines" && p.lut ~= "lines"
        hs.lut.UserData.lastLut = p.lut; % remember old lut
    end
    
    p.hsI(1).UserData.(cmd) = val;
    if any(cmd == ["lut" "lb" "ub" "volume"]), set_colorbar(hs); end
    if cmd == "volume", set_xyz(hs); end
    set_cdata(hs);
elseif cmd == "toggle" % turn on/off NIfTI
    i = evt.Indices(1);
    p = get_para(hs, i);
    if p.show == hs.files.Data{i,1}, return; end % no change
    p.show = ~p.show;
    p.hsI(1).UserData = p;
    
    states = {'off' 'on'};
    try %#ok<*TRYNC>
        set(p.hsI, 'Visible', states{p.show+1});
        if p.show, set_cdata(hs); end
        set_xyz(hs);
    end
elseif cmd == "mousemove"
    c = cell2mat(get(hs.ax(1:3), 'CurrentPoint'));
    c = c([1 3 5], 1:2); % 3x2
    x = cell2mat(get(hs.ax(1:3), 'XLim'));
    y = cell2mat(get(hs.ax(1:3), 'YLim'));
    I = cell2mat(get(hs.ijk, 'Value'))';
    if     c(1,1)>x(1,1) && c(1,1)<x(1,2) && c(1,2)>y(1,1) && c(1,2)<y(1,2) % sag
        I = [I(1) c(1,:)];
    elseif c(2,1)>x(2,1) && c(2,1)<x(2,2) && c(2,2)>y(2,1) && c(2,2)<y(2,2) % cor
        I = [c(2,1) I(2) c(2,2)];
    elseif c(3,1)>x(3,1) && c(3,1)<x(3,2) && c(3,2)>y(3,1) && c(3,2)<y(3,2) % tra
        I = [c(3,:) I(3)];
    end
    set_xyz(hs, I);
elseif cmd == "open" % open on current fig or new fig
    pName = hs.pref.UserData.openPath;
    [fname, pName] = uigetfile([pName '/*.nii; *.hdr;*.nii.gz; *.hdr.gz'], ...
        'Select a NIfTI to view', 'MultiSelect', 'on');
    if isnumeric(fname), return; end
    fname = strcat([pName '/'], fname);
    if strcmp(h.Label, 'Open in new window'), nii_viewer(fname);
    else, nii_viewer(fname, fh);
    end
    return;
elseif cmd == "add" % add overlay
    vars = evalin('base', 'who');
    is_nii = @(v)evalin('base', ...
        sprintf('isstruct(%s) && all(isfield(%s,{''hdr'',''img''}))', v, v));
    for i = numel(vars):-1:1, if ~is_nii(vars{i}), vars(i) = []; end; end
    if ~isempty(vars)
        a = listdlg('SelectionMode', 'single', 'ListString', vars, ...
            'ListSize', [300 100], 'CancelString', 'File dialog', ...
            'Name', 'Select a NIfTI in the list or click File dialog');
        if ~isempty(a), fname = evalin('base', vars{a}); end
    end
    
    pName = hs.pref.UserData.addPath;
    if strcmp(h.Label, 'Add aligned overlay')
        if ~exist('fname', 'var')
            [fname, pName] = uigetfile([pName '/*.nii; *.hdr;*.nii.gz;' ...
                '*.hdr.gz'], 'Select overlay NIfTI');
            if ~ischar(fname), return; end
            fname = fullfile(pName, fname);
        end
        [mtx, pName] = uigetfile([pName '/*.mat;*_warp.nii;*_warp.nii.gz'], ...
            'Select FSL mat file or warp file transforming the nii to background');
        if ~ischar(mtx), return; end
        mtx = fullfile(pName, mtx);
        addOverlay(fname, fh, mtx);
    else
        if ~exist('fname', 'var')
            [fname, pName] = uigetfile([pName '/*.nii; *.hdr;*.nii.gz;' ...
                '*.hdr.gz'], 'Select overlay NIfTI', 'MultiSelect', 'on');
            if ~ischar(fname) && ~iscell(fname), return; end
            fname = get_nii(strcat([pName filesep], fname));
        end
        addOverlay(fname, fh);
    end
    setpref('nii_viewer_para', 'addPath', pName);
elseif cmd == "closeAll" % close all overlays
    for j = size(hs.files.Data,1):-1:1
        p = get_para(hs, j);
        if p.hsI(1) == hs.hsI(1), continue; end
        delete(p.hsI); % remove image
        hs.files.Data(j,:) = [];
    end
    nii_viewer_cb(hs.files, struct('Indices', [1 2]), "file", fh);
    set_xyz(hs);
    table_width(hs.files);
elseif cmd == "close" % close selected overlay
    jf = -hs.stack.Value;
    p = get_para(hs, jf);
    if p.hsI(1) == hs.hsI(1), return; end % background
    delete(p.hsI); % 3 view
    hs.files.Data(jf,:) = [];
    jf = max(1, jf-1);
    nii_viewer_cb(hs.files, struct('Indices', [jf 2]), "file", fh);
    set_xyz(hs);
    table_width(hs.files);
elseif any(cmd == ["hdr" "ext" "essential"]) % show hdr, ext or essential
    jf = -hs.stack.Value;
    p = get_para(hs, jf);
    if cmd == "hdr"
        hdr = p.hdr0;
    elseif cmd == "ext"
        if ~isfield(p.nii, 'ext')
            uialert(fh, 'No extension for the selected NIfTI', 'Error');
            return;
        end
        hdr = {};
        for i = 1:numel(p.nii.ext)
            if ~isfield(p.nii.ext(i), 'edata_decoded'), continue; end
            hdr{end+1} = p.nii.ext(i).edata_decoded; %#ok
        end
        if isempty(hdr)
            uialert(fh, 'No known extension for the selected NIfTI', 'Error');
            return;
        elseif isscalar(hdr), hdr = hdr{1};
        end
    elseif cmd == "essential"
        hdr = nii_essential(p);
    end
    nam = hs.files.Data{jf, 2};
    if ~isstrprop(nam(1), 'alpha'), nam = ['x' nam]; end % like genvarname
    nam(~isstrprop(nam, 'alphanum')) = '_'; % make it valid for var name
    nam = [nam '_' cmd{1}];
    nam = strrep(nam, '__', '_');
    n = numel(nam); nm = namelengthmax;
    if n>nm, nam(nm-4:n-4) = ''; end
    assignin('base', nam, hdr);
    evalin('base', ['openvar ' nam]);
elseif cmd == "cross" % show/hide crosshairs and RAS labels
    if strcmp(h.Checked, 'on')
        h.Checked = 'off';
        set([hs.cross(:)' hs.ras hs.xyz], 'Visible', 'off');
    else
        h.Checked = 'on';
        set([hs.cross(:)' hs.ras hs.xyz], 'Visible', 'on');
    end
elseif cmd == "color" % crosshair color
    c = uisetcolor(hs.ras(1).Color, 'Pick crosshair color');
    if numel(c) ~= 3, return; end
    set([hs.cross(:)' hs.ras hs.xyz], 'Color', c);
    fh.Visible = 'off'; fh.Visible = 'on'; % bring to front?
elseif cmd == "thickness" % crosshair thickness
    set(hs.cross(:)', 'LineWidth', sscanf(h.Label, '%f'));
elseif cmd == "gap" % crosshair gap
    hs.gap = min(hs.pixdim) ./ hs.pixdim * sscanf(h.Label, '%f') / 2;
    guidata(fh, hs);
    set_cross(hs, 1:3);
elseif cmd == "copy" % copy figure into clipboard
    try res = hs.pref.UserData.dpi; catch, res = hs.hsN.pref.UserData.dpi; end
    if res==0, res = 150; end
    copygraphics(hs.frame, 'Resolution', res, 'BackgroundColor', 'current');
elseif cmd == "save" % save figure as picture
    [fname, pName] = uiputfile(['*.' h.Label], 'Input file name to save figure');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    try res = hs.pref.UserData.dpi; catch, res = hs.hsN.pref.UserData.dpi; end
    if res==0, res = 150; end
    exportgraphics(hs.frame, fname, 'Resolution', res, 'BackgroundColor', 'current');
elseif cmd == "colorbar" % colorbar on/off
    if strcmpi(hs.colorbar.Visible, 'on')
        hs.colorbar.Visible = 'off';
        h.Checked = 'off';
    else
        hs.colorbar.Visible = 'on';
        h.Checked = 'on';
        set_colorbar(hs);
    end
elseif cmd == "about"
    getVersion = dicm2nii('', 'getVersion', 'func_handle');
    str = sprintf(['nii_viewer.m by Xiangrui Li\n\n' ...
        'Last updated on %s\n\n', ...
        'Feedback to: xiangrui.li@gmail.com\n'], getVersion());
    helpdlg(str, 'About nii_viewer')
elseif cmd == "stack"
    st = hs.files.StyleConfigurations;
    jf = st.TargetIndex{1};
    n = size(hs.files.Data, 1);
    if h.Type == "uispinner"
        if -h.Value<jf, h.Tag = 'up'; else, h.Tag = 'down'; end
    end
    switch h.Tag
        case 'up' % one level up
            if jf==1, return; end
            ind = [1:jf-2 jf jf-1 jf+1:n]; jf = jf-1;
        case 'down' % one level down
            if jf==n, return; end
            ind = [1:jf-1 jf+1 jf jf+2:n]; jf = jf+1;
        case 'top'
            if jf==1, return; end
            ind = [jf 1:jf-1 jf+1:n]; jf = 1;
        case 'bottom'
            if jf==n, return; end
            ind = [1:jf-1 jf+1:n jf]; jf = n;
        otherwise
            error('Unknown stack level: %s', h.Tag);
    end
    
    im_stack(hs.ax, ind);
    hs.files.Data = hs.files.Data(ind,:);
    hs.stack.Value = -jf;
    hs.files.removeStyle; hs.files.addStyle(st.Style, 'row', jf);
    set_xyz(hs);
elseif cmd == "zoom"
    m = str2double(h.Label);
    a = min(hs.dim) / m;
    if a<1, m = min(hs.dim); end
    set_zoom(m, hs);
elseif cmd == "background"
    if strcmpi(h.Checked, 'on')
        h.Checked = 'off';
        hs.frame.BackgroundColor = [0 0 0];
        hs.colorbar.EdgeColor = [1 1 1];
    else
        h.Checked = 'on';
        hs.frame.BackgroundColor = [1 1 1];
        hs.colorbar.EdgeColor = [0 0 0];
    end
    set_cdata(hs);
elseif cmd == "flipLR"
    hs.pref.UserData.rightOnLeft = strcmpi(h.Checked, 'on');
    if hs.pref.UserData.rightOnLeft
        h.Checked = 'off';
        set(hs.ax([2 3]), 'XDir', 'normal');
        set(hs.ras([3 5]), 'String', 'L');
    else
        h.Checked = 'on';
        set(hs.ax([2 3]), 'XDir', 'reverse');
        set(hs.ras([3 5]), 'String', 'R');
    end
elseif cmd == "layout"
    layout = str2double(h.Tag);
    if h.Parent.UserData == layout, return; end
    h.Parent.UserData = layout;
    htP = hs.panel.Position(4);
    [siz, axPos, figPos] = plot_pos(hs.dim.*hs.pixdim, layout);
    fh.Position = [figPos siz+[2 htP]];
    fh.Visible = 'off'; drawnow; fh.Visible = 'on'; drawnow; % allow axPos update
    for i = 1:4, hs.ax(i).Position = axPos(i,:); end
elseif cmd == "keyHelp"
    str = sprintf([ ...
        'Key press available when focus is not in a number dialer:\n\n' ...
        'Left or Right arrow key: Move crosshair left or right.\n\n' ...
        'Up or Down arrow key: Move crosshair superior or inferior.\n\n' ...
        '[ or ] key: Move crosshair posterior or anterior.\n\n' ...
        '< or > key: Decrease or increase volume number.\n\n' ...
        'Ctrl + or - key: Zoom in or out by 10%% around crosshair.\n\n' ...
        'A: Toggle on/off crosshair.\n\n' ...
        'C: Crosshair to view center.\n\n' ...
        'Space: Toggle on/off selected image.\n\n' ...
        'F1: Show help text.\n']);
    helpdlg(str, 'Key Shortcut');
elseif cmd == "center" % image center
    p = get_para(hs);
    dim = p.nii.hdr.dim(2:4);
    c = round(hs.bg.Ri * (p.R * [dim/2-1 1]')) + 1;
    for i = 1:3, hs.ijk(i).Value = c(i); end
    nii_viewer_cb([], [], "ijk", fh);
elseif cmd == "viewCenter"
    c(3) = mean(hs.ax(1).YLim);
    c(1) = mean(hs.ax(2).XLim);
    c(2) = mean(hs.ax(1).XLim);
    c = double(round(c-0.5));
    for i = 1:3, hs.ijk(i).Value = c(i); end
    nii_viewer_cb([], [], "ijk", fh);
elseif cmd == "toXYZ"
    c0 = cell2mat(get(hs.ijk, 'Value'));
    c0 = hs.bg.R * [c0-1; 1];
    c0 = sprintf('%g %g %g', round(c0(1:3)));
    str = 'X Y Z coordinates in mm';
    while 1
        a = inputdlg(str, 'Crosshair to xyz', 1, {c0});
        if isempty(a), return; end
        c = sscanf(a{1}, '%g %g %g');
        if numel(c) == 3, break; end
    end
    c = round(hs.bg.Ri * [c(:); 1]) + 1;
    for i = 1:3, hs.ijk(i).Value = c(i); end
    nii_viewer_cb([], [], 'ijk', fh);
elseif cmd == "toValue"
    def = getappdata(h, 'Value');
    if isempty(def), def = 1; end
    def = num2str(def);
    str = 'Input the voxel value';
    while 1
        a = inputdlg(str, 'Crosshair to a value', 1, {def});
        if isempty(a), return; end
        val = sscanf(a{1}, '%g');
        if ~isnan(val), break; end
    end
    setappdata(h, 'Value', val);
    jf = -hs.stack.Value;
    p = get_para(hs, jf);
    img = p.nii.img(:,:,:, hs.volume.Value);
    c = find(img(:)==val, 1);
    if isempty(c)
        nam = strtok(hs.files.Data{jf,2}, '(');
        uialert(fh, sprintf('No value of %g found in %s', val, nam), 'Error');
        return;
    end
    dim = size(img); dim(numel(dim)+1:3) = 1;
    [c(1), c(2), c(3)] = ind2sub(dim, c); % ijk+1
    c = round(hs.bg.Ri * (p.R * [c(:)-1; 1])) + 1;
    for i = 1:3, hs.ijk(i).Value = c(i); end
    nii_viewer_cb([], [], 'ijk', fh);
elseif cmd == "cog" % crosshair to img COG
    p = get_para(hs);
    img = p.nii.img(:,:,:, hs.volume.Value);
    c = img_cog(img);
    if any(isnan(c)), uialert(fh, 'No valid COG found', 'Error'); return; end
    c = round(hs.bg.Ri * (p.R * [c-1; 1])) + 1;
    for i = 1:3, hs.ijk(i).Value = c(i); end
    nii_viewer_cb([], [], 'ijk', fh);
elseif cmd == "maximum" % crosshair to img max
    p = get_para(hs);
    img = p.nii.img(:,:,:,hs.volume.Value);
    if numel(unique(img(:))) < 2
        uialert(fh, 'All value are the same. No maximum!', 'Warning');
        return;
    end
    img = smooth23(img, 'gaussian', 5);
    img(isnan(img)) = 0;
    img = abs(img);
    [~, I] = max(img(:));
    dim = size(img); dim(end+1:3) = 1;
    c = zeros(3, 1);
    [c(1), c(2), c(3)] = ind2sub(dim, I);
    c = round(hs.bg.Ri * (p.R * [c-1; 1])) + 1;
    for i = 1:3, hs.ijk(i).Value = c(i); end
    nii_viewer_cb([], [], 'ijk', fh);
elseif cmd == "custom" % add custom lut
    p = get_para(hs);
    pName = fileparts(p.nii.hdr.file_name);
    [fname, pName] = uigetfile([pName '/*.lut'], 'Select LUT file for current overlay');
    if ~ischar(fname), return; end
    fid = fopen(fullfile(pName, fname));
    p.map = fread(fid, '*uint8');
    fclose(fid);
    if mod(numel(p.map),3)>0 || sum(p.map<8)<3 % guess text file
        try p.map = str2num(char(p.map'));
        catch, uialert(fh, 'Unrecognized LUT file', 'Error'); return;
        end
        if any(p.map(:)>1), p.map = single(p.map) / 255; end
    else
        p.map = reshape(p.map, [], 3);
        p.map = single(p.map) / 255;
    end
    if isequal(p.nii.img, round(p.nii.img))
        try p.map = p.map(1:max(p.nii.img(:))+1, :); end
    end
    p.lut = "custom";
    p.hsI(1).UserData = p;
    set(hs.lut, 'Value', p.lut, 'Enable', 'off');
    set_cdata(hs);
    set_colorbar(hs);
elseif cmd == "tc" % time course or std
    jf = -hs.stack.Value;
    p = get_para(hs, jf);
    nam = strtok(hs.files.Data{jf,2}, '(');
    
    labl = strrep(h.Label, ' ...', '');
    r = num2str(getappdata(h, 'radius'));
    r = inputdlg('Radius around crosshair (mm):', labl, 1, {r});
    if isempty(r), return; end
    r = str2double(r{1});
    setappdata(h, 'radius', r);
    c = [hs.ijk.Value]; % ijk for background
    c = hs.bg.R * [c-1 1]'; % in mm now
    c = c(1:3);
    
    b = xyzr2roi(c, r, p.nii.hdr); % overlay space
    img = p.nii.img;
    dim = size(img);
    img = reshape(img, [], prod(dim(4:end)))';
    img = img(:, b(:));
    if strcmp(labl, 'Time course')
        img = mean(single(img), 2);
    else
        img = std(single(img), [], 2);
    end
    fh1 = figure;
    plot(img);
    xlabel('Volume number');
    c = sprintf('(%g,%g,%g)', round(c));
    set(fh1, 'Name', [nam ' ' lower(labl) ' around voxel ' c]);
elseif cmd == "hist" % plot histgram
    jf = -hs.stack.Value;
    p = get_para(hs, jf);
    img = p.nii.img(:,:,:, hs.volume.Value);
    img = sort(img(:));
    img(isnan(img)) = [];
    img(img<hs.lb.Value) = [];
    img(img>hs.ub.Value) = [];
    nv = numel(img);
    img0 = unique(img);
    nu = numel(img0);
    n = max([nv/2000 nu/20 10]);
    n = min(round(n), nu);
    if n == nu, edges = img0;
    else, edges = linspace(0,1,n)*double(img(end)-img(1)) + double(img(1));
    end
    nam = strtok(hs.files.Data{jf,2}, '(');
    figure('NumberTitle', 'off', 'Name', nam);
    hh = histogram(img, edges, 'EdgeColor', "none");
    hh.BinCounts = hh.BinCounts / sum(hh.BinCounts) / median(diff(edges));
    xlabel('Voxel values'); ylabel('Probability density');
    title('Histogram between min and max values');
elseif cmd == "saveVolume" % save 1 or more volumes as a nifti
    p = get_para(hs);
    nam = p.nii.hdr.file_name;
    t = p.volume;
    while 1
        a = inputdlg('Volume indice to save (2:4 for example)', ...
            'Save Volume', 1, {num2str(t)});
        if isempty(a), return; end
        t = str2num(a{1});
        if ~isempty(t), break; end
    end
    pName = fileparts(nam);
    [fname, pName] = uiputfile([pName '/*.nii;*.nii.gz'], ...
        'Input name to save volume as');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    
    nii = nii_tool('load', nam); % re-load to be safe
    nii.img =  nii.img(:,:,:,t);
    nii_tool('save', nii, fname);
elseif cmd == "ROI" % save sphere
    c0 = cell2mat(get(hs.ijk, 'Value'));
    c0 = hs.bg.R * [c0-1; 1];
    c0 = sprintf('%g %g %g', round(c0(1:3)));
    str = {'X Y Z coordinates in mm' 'Radius in mm'};
    while 1
        a = inputdlg(str, 'Sphere ROI', 1, {c0 '6'});
        if isempty(a), return; end
        c = sscanf(a{1}, '%g %g %g');
        r = sscanf(a{2}, '%g');
        if numel(c) == 3, break; end
    end
    
    p = get_para(hs);
    pName = fileparts(p.nii.hdr.file_name);
    [fname, pName] = uiputfile([pName '/*.nii;*.nii.gz'], ...
        'Input file name to save ROI into');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    
    b = xyzr2roi(c, r, p.nii.hdr);
    p.nii.img = single(b); % single better supported by FSL
    nii_tool('save', p.nii, fname);
elseif cmd == "cropNeck"
    k0 = hs.ijk(3).Value - 1;
    p = get_para(hs);
    nam = p.nii.hdr.file_name;
    pName = fileparts(nam);
    [fname, pName] = uiputfile([pName '/*.nii;*.nii.gz'], ...
        'Input file name to save cropped image');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    
    R = nii_xform_mat(p.hdr0, hs.form_code); % original R
    k = hs.bg.Ri * R * grid3(p.hdr0.dim(2:4)); % background ijk
    
    nii = nii_tool('load', nam);
    d = size(nii.img);
    nii.img = reshape(nii.img, prod(d(1:3)), []);
    nii.img(k(3,:)<k0, :) = -nii.hdr.scl_inter / nii.hdr.scl_slope;
    nii.img = reshape(nii.img, d);
    nii_tool('save', nii, fname);
elseif cmd == "closeFig"
    try close(fh.UserData); end % cii_view
    delete(fh); return;
elseif cmd == "file"
    jf = evt.Indices(:,1);
    if numel(jf)>1 || evt.Indices(2)==1, return; end % multi-select or checkbox    
    p = get_para(hs, jf);
    nVol = size(p.nii.img, 4);
    
    a = hs.lut.UserData.lutStr; % show only meaningful lut for the file
    if p.lut ~= "custom", a(a == "custom") = []; end
    if nVol ~= 3 || max(abs(p.nii.img(:)))>1, a(a == "lines") = []; end
    if p.lut ~= "RGB" && nVol~=3 && size(p.nii.img,8)<3, a(a == "RGB") = []; end
    if ~startsWith(p.lut, "phase"), a(startsWith(a, "phase")) = []; end
    hs.lut.Items = a;
 
    st = h.StyleConfigurations;
    nf = size(h.Data,1);
    hs.stack.Enable = nf>1;
    if nf>1, hs.stack.Limits(1) = -nf; end
    hs.stack.Value = -jf;
    try h.Selection = [jf 2]; end % since R2020b
    h.removeStyle; h.addStyle(st.Style, 'row', jf);    
    cellfun(@(c)set(hs.(c), 'Value', p.(c)), {'lb' 'ub' 'alpha' 'volume' 'lut' 'smooth' 'interp'});
 
    hs.lb.Step = p.lb_step;
    hs.ub.Step = p.ub_step;
    set(hs.volume, 'Enable', nVol>1, 'ToolTip', "Volume number, 1:"+nVol);
    if nVol>1, hs.volume.Limits(2) = nVol; end
    set_colorbar(hs);
    
    off_on = {'off' 'on'};
    hs.interp.Enable = off_on{isfield(p, 'R0')+1};
    set(hs.overlay(1:4), 'Enable', off_on{(nf>1)+1}); % stack & Remove overlays
    hs.overlay(5).Enable = off_on{(p.hsI(1) ~= hs.hsI(1))+1}; % Remove overlay
    hs.lut.Enable = off_on{2-(p.lut=="custom" || size(p.nii.img,8)>1)};
        
elseif cmd == "drop"
    try
        nii = get_nii(evt.names);
        if evt.ctrlKey, addOverlay(nii, fh); % add overlay
        else, nii_viewer(nii, fh); return; % open
        end
    catch me
        uialert(fh, me.message, 'Error');
    end
elseif cmd == "MP4" % save slices as movie
    str = {'Which view to slice? 1:Sag; 2:Cor; 3:Tra', 'One view only?'};
    a = inputdlg(str, 'Export Movie', 1, {'3' 'No'});
    if isempty(a), return; end
    d = str2double(strtrim(a{1}));
    if ~any(d==1:3), errordlg('Input must be 1, 2 or 3'); return; end
    oneView = strcmpi('y', a{2}(1));
    str = {'Slice range (default all)' 'Movie frames per second'};
    a = inputdlg(str, 'Export Movie', 1, {['1:' num2str(hs.dim(d))] '4'});
    if isempty(a), return; end
    rg = eval(['[' a{1} ']']);
    fps = str2double(a{2});
    p = get_para(hs);
    pName = fileparts(p.nii.hdr.file_name);
    [fname, pName] = uiputfile([pName '/*.mp4'], 'Input file name to save the movie');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);

    if oneView, rect = round(getpixelposition(hs.ax(d)));
    else, rect = hs.frame.Position;
    end
    rect(3:4) = fix(rect(3:4)/2) * 2;
    vw = mp4_video(fname, fps, rect);
    for i = rg
        hs.ijk(d).Value = i;
        nii_viewer_cb(hs.ijk(d), [], "ijk", fh);
        vw.addFrame(fh);
    end
    vw.save();
else
    error('Unknown Callback: %s', cmd);
end
focus_off(fh);
try cii_view_cb(fh.UserData, [], cmd); end
 
%% zoom in/out with a factor
function set_zoom(m, hs)
c = hs.dim(:) / 2;
if m <= 1, I = c; % full view regardless of crosshair location
else, I = cell2mat(get(hs.ijk, 'Value'));
end
lim = round([I I] + c/m*[-1 1]) + 0.5;
axis(hs.ax(1), [lim(2,:) lim(3,:)]);
axis(hs.ax(2), [lim(1,:) lim(3,:)]);
axis(hs.ax(3), [lim(1,:) lim(2,:)]);
 
%% KeyPressFcn for figure
function KeyPressFcn(fh, evt)
if any(strcmp(evt.Key, evt.Modifier)), return; end % only modifier
hs = guidata(fh);
if ~isempty(intersect({'control' 'command'}, evt.Modifier))
    switch evt.Key
        case {'add' 'equal'}
            [dim, i] = min(hs.dim);
            if     i==1, d = hs.ax(2).XLim;
            elseif i==2, d = hs.ax(1).XLim;
            else,        d = hs.ax(1).YLim;
            end
            d = abs(diff(d'));
            if d<=3, return; end % ignore
            m = dim / d * 1.1;
            if round(dim/2/m)==d/2, m = dim / (d-1); end
            set_zoom(m, hs);
        case {'subtract' 'hyphen'}
            d = abs(diff(hs.ax(2).XLim));
            m = hs.dim(1) / d;
            if m<=1, return; end
            m = m / 1.1;
            if round(hs.dim(1)/2/m)==d/2, m = hs.dim(1) / (d+1); end
            if m<1.01, m = 1; end
            set_zoom(m, hs);
    end
    return;
end
 
switch evt.Key
    case 'leftarrow'
        if hs.ijk(1).Value==1, return; end
        hs.ijk(1).Value = hs.ijk(1).Value - 1;
        nii_viewer_cb(hs.ijk(1), [], 'ijk', fh)
    case 'rightarrow'
        if hs.ijk(1).Value==hs.dim(1), return; end
        hs.ijk(1).Value = hs.ijk(1).Value + 1;
        nii_viewer_cb(hs.ijk(1), [], 'ijk', fh)
    case 'uparrow'
        if hs.ijk(3).Value==hs.dim(3), return; end
        hs.ijk(3).Value = hs.ijk(3).Value + 1;
        nii_viewer_cb(hs.ijk(3), [], 'ijk', fh)
    case 'downarrow'
        if hs.ijk(3).Value==1, return; end
        hs.ijk(3).Value = hs.ijk(3).Value - 1;
        nii_viewer_cb(hs.ijk(3), [], 'ijk', fh)
    case 'rightbracket' % ]
        if hs.ijk(2).Value==hs.dim(2), return; end
        hs.ijk(2).Value = hs.ijk(2).Value + 1;
        nii_viewer_cb(hs.ijk(2), [], 'ijk', fh)
    case 'leftbracket' % [
        if hs.ijk(2).Value==1, return; end
        hs.ijk(2).Value = hs.ijk(2).Value - 1;
        nii_viewer_cb(hs.ijk(2), [], 'ijk', fh)
    case 'period' % . or >
        if hs.volume.Value==round(hs.volume.Limits(2)), return; end
        hs.volume.Value = hs.volume.Value + 1;
        nii_viewer_cb(hs.volume, [], 'volume', fh)
    case 'comma' % , or <
        if hs.volume.Value==1, return; end
        hs.volume.Value = hs.volume.Value - 1;
        nii_viewer_cb(hs.volume, [], 'volume', fh)
    case 'c'
        nii_viewer_cb([], [], 'viewCenter', hs.fig);
    case {'x' 'space'}
        i = -hs.stack.Value;
        hs.files.Data{i,1} = ~hs.files.Data{i,1};
        nii_viewer_cb(hs.files, struct('Indices', [i 1]), 'toggle', hs.fig);
    case 'a'
        h = findobj(hs.fig, 'Type', 'uimenu', 'Label', 'Show crosshair');
        nii_viewer_cb(h, [], 'cross', hs.fig);
    case 'f1'
        doc nii_viewer;
    case 'tab' % prevent tab from cycling ui components
        focus_off(fh);
end
 
%% update CData/AlphaData for 1 or 3 of the sag/cor/tra views
function set_cdata(hs, iAxis)
if nargin<2, iAxis = 1:3; end
nFile = size(hs.files.Data, 1);
for i = 1:nFile
    p = get_para(hs, i);
    if ~p.show, continue; end % save time, but need to update when enabled
    if p.lut == "lines" % special case: do it separately
        vector_lines(hs, i, iAxis); continue; 
    elseif ~strcmpi(p.hsI(1).Type, 'image') % was "lines"
        delete(p.hsI); % delete quiver
        p.hsI = copyimg(hs);
        p.hsI(1).UserData = p; % update whole UserData
        if i>1, im_stack(hs.ax, [2:i 1 i+1:nFile]); end
    end
    t = round(p.volume);
    img = p.nii.img;
    isRGB = size(img, 8)>2;
    if isRGB % avoid indexing for single vol img: could speed up a lot
        img = permute(img(:,:,:,t,:,:,:,:), [1:3 8 4:7]);
    elseif size(img,4)>1 && p.lut~="RGB"
        img = img(:,:,:,t);
    end    
    if ~isfloat(img)
        img = single(img);
        if isfield(p, 'scl_slope')
            img = img * p.scl_slope + p.scl_inter;
        end
    end
    
    if isfield(p, 'mask'), img = img .* p.mask; end
    if isfield(p, 'modulation'), img = img .* p.modulation; end
 
    if startsWith(p.lut, "phase") % interp/smooth both mag and phase
        img(:,:,:,2) = p.phase(:,:,:,t);
    end
    
    dim4 = size(img, 4);
    for ix = iAxis
        ind = round(hs.ijk(ix).Value);
        if ind<1 || ind>hs.dim(ix), continue; end
        ii = {':' ':' ':'};
        io = ii;
        d = hs.dim;
        d(ix) = 1; % 1 slice at dim ix
        im = zeros([d dim4], 'single');
        
        if isfield(p, 'R0') % interp, maybe smooth
            I = grid3(d);
            I(ix,:) = ind-1;
            
            if isfield(p, 'warp')
                iw = {':' ':' ':' ':'}; iw{ix} = ind;
                warp = p.warp(iw{:});
                warp = reshape(warp, [], 3)'; warp(4,:) = 0;
                I = p.Ri * (p.R0 * I + warp) + 1;
            else
                I = p.Ri * p.R0 * I + 1; % ijk+1 for overlay img
            end
            
            for j = 1:dim4
                if p.smooth
                    ns = 3; % number of voxels (odd) used for smooth
                    d3 = d; d3(ix) = ns; % e.g. 3 slices
                    b = zeros(d3, 'single');
                    I0 = I(1:3,:);
                    for k = 1:ns % interp for each slice
                        I0(ix,:) = I(ix,:) - (ns+1)/2 + k;
                        a = interp3a(img(:,:,:,j), I0, p.interp);
                        ii{ix} = k; b(ii{:}) = reshape(a, d);
                    end
                    b = smooth23(b, 'gaussian', ns);
                    io{ix} = (ns+1)/2;
                    im(:,:,:,j) = b(io{:}); % middle one
                else
                    a = interp3a(img(:,:,:,j), I, p.interp);
                    im(:,:,:,j) = reshape(a, d);
                end
            end
        elseif p.smooth % smooth only
            ns = 3; % odd number of slices to smooth
            ind1 = ind - (ns+1)/2 + (1:ns); % like ind+(-1:1)
            if any(ind1<1 | ind1>hs.dim(ix)), ind1 = ind; end % 2D
            ii{ix} = ind1;
            io{ix} = mean(1:numel(ind1)); % middle slice
            for j = 1:dim4
                a = smooth23(img(ii{:},j), 'gaussian', ns);
                im(:,:,:,j) = a(io{:});
            end
        else % no interp or smooth
            io{ix} = ind;
            im(:) = img(io{:}, :);
        end
        
        if     ix == 1, im = permute(im, [3 2 4 1]);
        elseif ix == 2, im = permute(im, [3 1 4 2]);
        elseif ix == 3, im = permute(im, [2 1 4 3]);
        end
        
        if ~isRGB % not NIfTI RGB
            [im, alfa] = lut2img(im, p);
        elseif dim4 == 3 % NIfTI RGB
            if max(im(:))>2, im = im / 255; end % guess uint8
            im(im>1) = 1; im(im<0) = 0;
            alfa = sum(im,3) / dim4; % avoid mean
        elseif dim4 == 4 % NIfTI RGBA
            if max(im(:))>2, im = im / 255; end % guess uint8
            im(im>1) = 1; im(im<0) = 0;
            alfa = im(:,:,4);
            im = im(:,:,1:3);
        else
            error('Unknown data type: %s', p.nii.hdr.file_name);
        end
        
        if p.hsI(1) == hs.hsI(1) && isequal(hs.frame.BackgroundColor, [1 1 1])
            alfa = img2mask(alfa);
        elseif dim4 ~= 4
            alfa = alfa > 0;
        end
        alfa = p.alpha * alfa;
        set(p.hsI(ix), 'CData', im, 'AlphaData', alfa);
    end
end
 
%% Add an overlay
function addOverlay(fname, fh, mtx)
hs = guidata(fh);
frm = hs.form_code;
aligned = nargin>2;
R_back = hs.bg.R;
R0 = nii_xform_mat(hs.bg.hdr, frm(1)); % original background R
[~, perm, flp] = reorient(R0, hs.bg.hdr.dim(2:4), 0);
if aligned % aligned mtx: do it in special way
    [p, ~, rg, dim] = read_nii(fname, frm, 0); % no re-orient
    try
        if any(regexpi(mtx, '\.mat$'))
            R = load(mtx, '-ascii');
            if ~isequal(size(R), [4 4])
                error('Invalid transformation matrix file: %s', mtx);
            end
        else % see nii_xform
            R = eye(4);
            warp = nii_tool('img', mtx); % FSL warp nifti
            if ~isequal(size(warp), [hs.bg.hdr.dim(2:4) 3])
                error('warp file and template file img size don''t match.');
            end
            if det(R0(1:3,1:3))<0, warp(:,:,:,1) = -warp(:,:,:,1); end
            if ~isequal(perm, 1:3)
                warp = permute(warp, [perm 4]);
            end
            for j = 1:3
                if flp(j), warp = flip(warp, j); end
            end
            p.warp = warp;
            p.R0 = R_back; % always interp
        end
    catch me
        uialert(fh, me.message, 'Error');
        return;
    end
 
    % see nii_xform for more comment on following method
    R = R0 / diag([hs.bg.hdr.pixdim(2:4) 1]) * R * diag([p.pixdim 1]);
    [~, i1] = max(abs(p.R(1:3,1:3)));
    [~, i0] = max(abs(R(1:3,1:3)));
    flp = sign(R(i0+[0 4 8])) ~= sign(p.R(i1+[0 4 8]));
    if any(flp)
        rotM = diag([1-flp*2 1]);
        rotM(1:3,4) = (dim-1).* flp;
        R = R / rotM;
    end
            
    [p.R, perm, p.flip] = reorient(R, dim); % in case we apply mask to it
    if ~isequal(perm, 1:3)
        dim = dim(perm);
        p.pixdim = p.pixdim(perm);
        p.nii.img = permute(p.nii.img, [perm 4:8]);
    end
    for j = 1:3
        if p.flip(j), p.nii.img = flip(p.nii.img, j); end
    end
    p.alignMtx = mtx; % info only for NIfTI essentials
else % regular overlay
    [p, frm, rg, dim] = read_nii(fname, frm);
 
    % Silently use another background R_back matching overlay: very rare
    if frm>0 && frm ~= hs.form_code(1) && any(frm == hs.form_code)
        R = nii_xform_mat(hs.bg.hdr, frm); % img alreay re-oriented
        R_back = reorient(R, hs.bg.hdr.dim(2:4));
    elseif frm==0 && isequal(p.hdr0.dim(2:4), hs.bg.hdr.dim(2:4))
        p.R = hs.bg.R; % best guess: use background xform
        p.perm = perm;
        p.nii.img = permute(p.nii.img, [p.perm 4:8]);
        for i = 1:3
            if p.flip(i) ~= flp(i)
                p.nii.img = flip(p.nii.img, i);
            end
        end
        p.flip = flp;
        uialert(fh, ['There is no valid coordinate system for the overlay. ' ...
         'The viewer supposes the same coordinate as the background.'], ...
         'Missing valid tranformation');
    elseif frm ~= 2 && ~any(frm == hs.form_code)
        uialert(fh, ['There is no matching coordinate system between the overlay ' ...
         'image and the background image. The overlay is likely meaningless.'], ...
         'Transform Inconsistent');
    end
end
 
singleVol = 0;
nv = size(p.nii.img, 4);
if nv>1 && numel(p.nii.img)>1e7 % load all or a single volume
    if isfield(p.nii, 'NamedMap')
        nams = cell(1, numel(p.nii.NamedMap));
        for i = 1:numel(nams), nams{i} = p.nii.NamedMap{i}.MapName; end
        a = listdlg('PromptString', 'Load "All" or one of the map:', ...
            'SelectionMode', 'single', 'ListString', ['All' nams]);
        if a==1, a = 'All'; else, a = a - 1; end
    else
        str = ['Input ''all'' or a number from 1 to ' num2str(nv)];
        while 1
            a = inputdlg(str, 'Volumes to load', 1, {'all'});
            if isempty(a), return; end
            a = strtrim(a{1});
            if ~isstrprop(a, 'digit'), break; end
            a = str2num(a);
            if isequal(a,1:nv) || (isscalar(a) && a>=1 && a<=nv && mod(a,1)==0)
                break;
            end
        end
    end
    if isnumeric(a) && isscalar(a)
        singleVol = a;
        p.nii.img = p.nii.img(:,:,:,a);
        if isfield(p.nii, 'cii')
            p.nii.cii{1} = p.nii.cii{1}(:,a);
            p.nii.cii{2} = p.nii.cii{2}(:,a);
        end
        if isfield(p.nii, 'NamedMap')
            try    p.nii.NamedMap = p.nii.NamedMap(a);
            catch, p.nii.NamedMap = p.nii.NamedMap(1);
            end
            try p.map = p.nii.NamedMap{1}.map; end
        end
        rg = get_range(p.nii);
    end
end
 
ii = [1 6 11 13:15]; % diag and offset: avoid large ratio due to small value
if ~isequal(hs.dim, dim) || any(abs(R_back(ii)./p.R(ii)-1) > 0.01)
    p.R0 = R_back;
end
p.Ri = inv(p.R);
 
if ~isreal(p.nii.img)
    p.phase = angle(p.nii.img); % -pi to pi
    p.phase = mod(p.phase/(2*pi), 1); % 0~1
    p.nii.img = abs(p.nii.img); % real now
end
 
hsI = findobj(hs.ax(1), 'Type', 'image', '-or', 'Type', 'quiver');
luts = string(arrayfun(@(a)a.UserData.lut, hsI, 'UniformOutput', false));
 
p.hsI = copyimg(hs); % duplicate image obj for overlay: will be at top
p.lb = rg(1); p.ub = rg(2);
p = dispPara(p, luts);
p.hsI(1).UserData = p;
 
[pName, niiName, ext] = fileparts(p.nii.hdr.file_name);
if strcmpi(ext, '.gz'), [~, niiName] = fileparts(niiName); end
if aligned, niiName = [niiName '(aligned)']; end
if singleVol, niiName = [niiName '(' num2str(singleVol) ')']; end
 
hs.files.Data = [{true niiName}; hs.files.Data];
st = hs.files.StyleConfigurations;
hs.files.removeStyle; hs.files.addStyle(st.Style, 'row', st.TargetIndex{1}+1);
nii_viewer_cb(hs.files, struct('Indices', [1 2]), "file", fh);
hs.pref.UserData.addPath = pName;
table_width(hs.files);
 
set_cdata(hs);
set_xyz(hs);
figure(hs.fig);
if isfield(p.nii, 'cii'), cii_view(hs); end
 
%% Reorient 4x4 R
function [R, perm, flp] = reorient(R, dim, leftHand)
% [R, perm, flip] = reorient(R, dim, leftHand)
% Re-orient transformation matrix R (4x4), so it will be diagonal major and
% positive at diagonal, unless the optional third input is true, which requires
% left-handed matrix, where R(1,1) will be negative. 
% The second input is the img space dimension (1x3). 
% The perm output, like [1 2 3] or a permutation of it, indicates if input R was
% permuted for 3 axis. The third output, flip (1x3 logical), indicates an axis 
% (AFTER perm) is flipped if true.
a = abs(R(1:3,1:3));
[~, ixyz] = max(a);
if ixyz(2) == ixyz(1), a(ixyz(2),2) = -1; [~, ixyz(2)] = max(a(:,2)); end
if any(ixyz(3) == ixyz(1:2)), ixyz(3) = setdiff(1:3, ixyz(1:2)); end
[~, perm] = sort(ixyz);
R(:,1:3) = R(:,perm);
flp = diag(R(1:3, 1:3))' < 0;
if nargin>2 && leftHand, flp(1) = ~flp(1); end
rotM = diag([1-flp*2 1]);
rotM(1:3, 4) = (dim(perm)-1) .* flp; % 0 or dim-1
R = R / rotM; % xform matrix after flip
 
%% Load, re-orient nii, extract essential nii stuff
% nifti may be re-oriented, p.hdr0 stores original nii.hdr
function [p, frm, rg, dim] = read_nii(fname, ask_code, reOri)
if nargin<2, ask_code = []; end
if ischar(fname), p.nii = nii_tool('load', fname);
else, p.nii = fname; fname = p.nii.hdr.file_name;
end
p.hdr0 = p.nii.hdr; % original hdr
c = p.nii.hdr.intent_code;
if c>=3000 && c<=3099 && isfield(p.nii, 'ext') && any([p.nii.ext.ecode] == 32)
    p.nii = cii2nii(p.nii);
end
 
if nargin<3 || reOri
    [p.nii, p.perm, p.flip] = nii_reorient(p.nii, 0, ask_code);
else
    p.perm = 1:3;
    p.flip = false(1,3);
end
 
dim = p.nii.hdr.dim(2:8);
dim(dim<1 | mod(dim,1)~=0) = 1;
if p.nii.hdr.dim(1)>4 % 4+ dim, put all into dim4
    if sum(dim(4:7)>1)>1
        uialert(fh, [fname ' has 5 or more dimension. Dimension above 4 are ' ...
            'all treated as volumes for visualization'], 'Warning');        
    end
    dim(4) = prod(dim(4:7)); dim(5:7) = 1;
    p.nii.img = reshape(p.nii.img, [dim size(p.nii.img, 8)]);
end
 
[p.R, frm] = nii_xform_mat(p.nii.hdr, ask_code);
dim = dim(1:3);
p.pixdim = p.nii.hdr.pixdim(2:4);
 
if size(p.nii.img,4)<4 && ~isfloat(p.nii.img)
    p.nii.img = single(p.nii.img); 
end
if p.nii.hdr.scl_slope==0, p.nii.hdr.scl_slope = 1; end
if p.nii.hdr.scl_slope~=1 || p.nii.hdr.scl_inter~=0
    if isfloat(p.nii.img)
        p.nii.img = p.nii.img * p.nii.hdr.scl_slope + p.nii.hdr.scl_inter;
    else
        p.scl_slope = p.nii.hdr.scl_slope;
        p.scl_inter = p.nii.hdr.scl_inter;
    end
end
 
% check if ROI labels available: the same file name with .txt extension
if c == 1002 % Label
    [pth, nam, ext] = fileparts(p.nii.hdr.file_name);
    nam1 = fullfile(pth, [nam '.txt']);
    if strcmpi(ext, '.gz') && ~exist(nam1, 'file')
        [~, nam] = fileparts(nam);
        nam1 = fullfile(pth, [nam '.txt']);
    end
    if exist(nam1, 'file') % each line format: 1 ROI_1
        fid = fopen(nam1);
        while 1
            ln = fgetl(fid);
            if ~ischar(ln), break; end
            [ind, a] = strtok(ln);
            ind = str2double(ind);
            try p.labels{ind} = strtrim(a); catch, end
        end
        fclose(fid);
    end
end
rg = get_range(p.nii, isfield(p, 'labels'));
try p.map = p.nii.NamedMap{1}.map; end
 
%% Return xform mat and form_code: form_code may have two if not to ask_code
function [R, frm] = nii_xform_mat(hdr, ask_code)
% [R, form] = nii_xform_mat(hdr, asked_code);
% Return the transformation matrix from a NIfTI hdr. By default, this returns
% the sform if available. If the optional second input, required form code, is
% provided, this will try to return matrix for that form code. The second
% optional output is the form code of the actually returned matrix.
fs = [hdr.sform_code hdr.qform_code]; % sform preferred
if fs(1)==fs(2), fs = fs(1); end % sform if both are the same
f = fs(fs>=1 & fs<=4); % 1/2/3/4 only
if isempty(f) || ~strncmp(hdr.magic, 'n', 1) % treat it as Analyze
    frm = 0;
    try % try spm style Analyze
        [pth, nam, ext] = fileparts(hdr.file_name);
        if strcmpi(ext, '.gz'), [~, nam] = fileparts(nam); end
        R = load(fullfile(pth, [nam '.mat']));
        R = R.M;
    catch % make up R for Analyze: suppose xyz order with left storage 
        R = [diag(hdr.pixdim(2:4)) -(hdr.dim(2:4).*hdr.pixdim(2:4)/2)'; 0 0 0 1];
        R(1,:) = -R(1,:); % use left handed
    end
    return;
end
 
if isscalar(f) || nargin<2 || isempty(ask_code) % only 1 avail or no ask_code
    frm = f;
else % numel(f) is 2, numel(ask_code) can be 1 or 2
    frm = f(f == ask_code(1));
    if isempty(frm) && numel(ask_code)>1, frm = f(f == ask_code(2)); end
    if isempty(frm) && any(f==2), frm = 2; end % use confusing code 2
    if isempty(frm), frm = f(1); end % no match to ask_code, use sform
end
 
if frm(1) == fs(1), R = [hdr.sform_mat; 0 0 0 1]; % match sform_code or no match
else, R = [hdr.qform_mat; 0 0 0 1]; % match qform_code
end
 
%% Estimate lower and upper bound of img display
function rg = get_range(nii, isLabel)
if size(nii.img, 8)>2 || any(nii.hdr.datatype == [128 511 2304]) % RGB / RGBA
    if max(nii.img(:))>2, rg = [0 255]; else, rg = [0 1]; end
    return;
elseif nii.hdr.cal_max~=0 && nii.hdr.cal_max>min(nii.img(:))
    rg = [nii.hdr.cal_min nii.hdr.cal_max];
    return;
end
 
img = nii.img(:,:,:,1);
img = img(:);
img(isnan(img) | isinf(img)) = [];
if ~isreal(img), img = abs(img); end
if ~isfloat(img)
    slope = nii.hdr.scl_slope; if slope==0, slope = 1; end
    img = single(img) * slope + nii.hdr.scl_inter;
end
img = double(img);
 
mi = min(img); ma = max(img);
if nii.hdr.intent_code > 1000 || (nargin>1 && isLabel)
    rg = [mi ma]; return;
elseif nii.hdr.intent_code == 2 % correlation
    rg = [0.3 1]; return;
end
 
ind = abs(img)>50;
if sum(ind)<numel(img)/10, ind = abs(img)>std(img)/2; end
im = img(ind);
mu = mean(im);
sd = std(im);
rg = mu + [-2 2]*sd;
if rg(1)<=0, rg(1) = sd/5; end
if rg(1)<mi || isnan(rg(1)), rg(1) = mi; end
if rg(2)>ma || isnan(rg(2)), rg(2) = ma; end
if rg(1)==rg(2), rg(1) = mi; if rg(1)==rg(2), rg(1) = 0; end; end
rg = round(rg, 2, 'significant'); % since 2014b
if rg(1)==rg(2), rg(1) = mi; end
if abs(rg(1))>10, rg(1) = floor(rg(1)/2)*2; end % even number
if abs(rg(2))>10, rg(2) = ceil(rg(2)/2)*2; end % even number
if mi<0 && abs(mi)>ma/2, rg(1) = -rg(2); end % like phasediff 
 
%% Draw vector lines, called by set_cdata
function vector_lines(hs, i, iaxis)
p = get_para(hs, i);
d = single(size(p.nii.img));
pixdim = hs.bg.hdr.pixdim(2:4); % before reorient
if strcmpi(p.hsI(1).Type, 'image') % just switched to "lines"
    cb = hs.hsI(1).ButtonDownFcn;
    delete(p.hsI);
    lut = hs.lut.UserData.lastLut;
    clr = lut2map(p, lut); clr = clr(end,:);
    for j = 1:3
        p.hsI(j) = quiver(hs.ax(j), 1, 1, 0, 0, 'Color', clr, ...
            'ShowArrowHead', 'off', 'AutoScale', 'off', 'ButtonDownFcn', cb);
    end
    crossFront(hs); % to be safe before next
    if i>1, im_stack(hs.ax, [2:i 1 i+1:size(hs.files.Data,1)]); end
    
    if isfield(p, 'R0') && ~isfield(p, 'ivec')
        I = p.R0 \ (p.R * grid3(d)) + 1;
        p.ivec = reshape(I(1:3,:)', d);
 
        R0 = normc(p.R0(1:3, 1:3));
        R  = normc(p.R(1:3, 1:3));
        [pd, j] = min(p.pixdim);
        p.Rvec = R0 / R * pd / pixdim(j);
    end
    p.hsI(1).UserData = p;
end
 
img = p.nii.img;
% This is needed since vec is in image ref, at least for fsl
img(:,:,:,p.flip) = -img(:,:,:,p.flip);
if isfield(p, 'mask') % ignore modulation
    img = img .* p.mask;
end
if any(abs(diff(pixdim))>1e-4) % non isovoxel background
    pd = pixdim / min(pixdim);
    for j = 1:3, img(:,:,:,j) = img(:,:,:,j) / pd(j); end
end
 
if isfield(p, 'Rvec')
    img = reshape(img, [], d(4));
    img = img * p.Rvec;
    img = reshape(img, d);
end
 
for ix = iaxis
    I = round(hs.ijk(ix).Value);
    j = 1:3; j(ix) = [];
    if isfield(p, 'ivec')
        I = abs(p.ivec(:,:,:,ix) - I);
        [~, I] = min(I, [], ix);
        
        ii = {1:d(1) 1:d(2) 1:d(3)};
        ii{ix} = single(1);
        [ii{1}, ii{2}, ii{3}] = ndgrid(ii{:});
        ii{ix} = single(I);
        io = {':' ':' ':' ':'}; io{ix} = 1;
        im = img(io{:});
        for k = 1:2
            im(:,:,:,k) = interp3(img(:,:,:,j(k)), ii{[2 1 3]}, 'nearest');
        end
    
        ind = sub2ind(d(1:3), ii{:});
        X = p.ivec(:,:,:,j(1)); X = permute(X(ind), [j([2 1]) ix]);
        Y = p.ivec(:,:,:,j(2)); Y = permute(Y(ind), [j([2 1]) ix]);    
    else
        ii = {':' ':' ':'};
        ii{ix} = I;
        im = img(ii{:}, j);
        [Y, X] = ndgrid(1:d(j(2)), 1:d(j(1)));
    end
    
    im = permute(im, [j([2 1]) 4 ix]);
    im(im==0) = nan; % avoid dots in emf and eps
    X = X - im(:,:,1)/2;
    Y = Y - im(:,:,2)/2;
    set(p.hsI(ix), 'XData', X, 'YData', Y, 'UData', im(:,:,1), 'VData', im(:,:,2));
end
 
%% Bring cross and label to front
function crossFront(hs)
for i = 1:3
    ch = hs.ax(i).Children;
    typ = arrayfun(@(a)a.Type, ch, 'UniformOutput', false);
    ind = strcmp(typ, 'line') | strcmp(typ, 'text');
    hs.ax(i).Children = [ch(ind); ch(~ind)];
end
 
%% Compute color map for LUT
function map = lut2map(p, lut)
if isfield(p, 'map')
    if isfield(p.nii, 'NamedMap')
        try map = p.nii.NamedMap{p.volume}.map; end
    else, map = p.map;
    end
    return;
end
 
if nargin<2, lut = p.lut; end
map = linspace(0,1,64)'*[1 1 1]; % gray
if     lut == "grayscale", return;
elseif lut == "red",    map(:,2:3) = 0;
elseif lut == "green",  map(:,[1 3]) = 0;
elseif lut == "blue",   map(:,1:2) = 0;
elseif lut == "violet", map(:,2) = 0;
elseif lut == "yellow", map(:,3) = 0;
elseif lut == "cyan",   map(:,1) = 0;
elseif any(lut == ["red-yellow" "autumn" "phase"]), map(:,3) = 0; map(:,1) = 1;
elseif lut == "blue-green", map(:,1) = 0; map(:,3) = flip(map(:,3));
elseif lut == "two-sided"
    map = map(1:2:end,:); % half
    map_neg = map;
    map(:,3) = 0; map(:,1) = 1; % red_yellow
    map_neg(:,1) = 0; map_neg(:,3) = flip(map_neg(:,3)); % blue_green
    map = [flip(map_neg,1); map];
elseif lut == "lines", map(:,2:3) = 0;
elseif lut == "phase3" % red-yellow-green-yellow-red
    a = map(:,1);
    map(1:32,3) = 0; map(1:16,1) = 1; map(17:32,2) = 1;
    map(1:16,2) = a(1:4:64); map(17:32,1) = a(64:-4:1);
    map(33:64,:) = map(32:-1:1,:);
elseif lut == "phase6" % red-yellow-green/violet-blue-cyan
    a = map(:,1);
    map(1:32,3) = 0; map(1:16,1) = 1; map(17:32,2) = 1;
    map(1:16,2) = a(1:4:64); map(17:32,1) = a(64:-4:1);
    map(33:48,2) = 0; map(33:48,3) = 1; map(33:48,1) = a(64:-4:1);
    map(49:64,1) = 0; map(49:64,3) = 1; map(49:64,2) = a(1:4:64);
elseif lut == "RGB"
else % matlab LUT
    map = feval(lut, 64);
end
 
%% Preference dialog
function pref_dialog(h, ~)
pf = getpref('nii_viewer_para');
d = dialog('Name', 'Preferences', 'Visible', 'off');
pos = getpixelposition(d);
pos(3:4) = [396 332];
hs.fig = ancestor(h, 'figure');
 
uicontrol(d, 'Style', 'text', 'Position', [8 306 300 22], ...
    'String', 'Default background image folder:', 'HorizontalAlignment', 'left');
hs.openPath = uicontrol(d, 'Style', 'edit', 'String', pf.openPath, ...
    'Position', [8 288 350 22], 'BackgroundColor', 'w', 'HorizontalAlignment', 'left', ...
    'TooltipString', 'nii_viewer will point to this folder when you "Open" image');
uicontrol(d, 'Position', [358 288 30 22], 'Tag', 'browse', ...
    'String', '...', 'Callback', @pref_dialog_cb);
 
hs.rightOnLeft = uicontrol(d, 'Style', 'popup', 'BackgroundColor', 'w', ...
    'Position', [8 252 380 22], 'Value', pf.rightOnLeft+1, ...
    'String', {'Neurological orientation (left on left side)' ...
               'Radiological orientation (right on left side)'}, ...
    'TooltipString', 'Display convention also applies to future use');
 
uicontrol(d, 'Style', 'text', 'Position', [8 210 40 22], ...
    'String', 'Layout', 'HorizontalAlignment', 'left', ...
    'TooltipString', 'Layout for three views');
 
% iconsFolder = fullfile(matlabroot,'/toolbox/matlab/icons/');
% iconUrl = strrep(['file:/' iconsFolder 'matlabicon.gif'],'\','/');
% str = ['<html><img src="' iconUrl '"/></html>'];
hs.layout = uicontrol(d, 'Style', 'popup', 'BackgroundColor', 'w', ...
    'Position', [50 214 338 22], 'Value', pf.layout, ...
    'String', {'one-row' 'two-row sag on right' 'two-row sag on left'}, ...
    'TooltipString', 'Layout for three views');
 
hs.mouseOver = uicontrol(d, 'Style', 'checkbox', ...
    'Position', [8 182 380 22], 'Value', pf.mouseOver, ...
    'String', 'Show coordinates and intensity when mouse moves over image', ...
    'TooltipString', 'Also apply to future use');
 
uipanel(d, 'Units', 'Pixels', 'Position', [4 110 390 56], ...
    'Title', 'For "Save NIfTI as" if interpolation is applicable');
str = {'nearest' 'linear' 'cubic' 'spline'};
val = find(strcmp(str, pf.interp));
uicontrol(d, 'Style', 'text', 'Position', [8 116 140 22], ...
    'String', 'Interpolation method:', 'HorizontalAlignment', 'right');
hs.interp = uicontrol(d, 'Style', 'popup', 'String', str, ...
    'Position', [150 120 68 22], 'Value', val, 'BackgroundColor', 'w');
 
uicontrol(d, 'Style', 'text', 'Position', [230 116 90 22], ...
    'String', 'Missing value:', 'HorizontalAlignment', 'right');
hs.extraV = uicontrol(d, 'Style', 'edit', 'String', num2str(pf.extraV), ...
    'Position', [324 120 60 22], 'BackgroundColor', 'w', ...
    'TooltipString', 'NaN or 0 is typical, but can be any number');
 
str = strtrim(cellstr(num2str([0 120 150 200 300 600 1200]')));
val = find(strcmp(str, pf.dpi));
uipanel(d, 'Units', 'Pixels', 'Position', [4 40 390 56], ...
    'Title', 'For "Save figure as" and "Copy figure"');
uicontrol(d, 'Style', 'text', 'Position', [8 46 90 22], ...
    'String', 'Resolution:', 'HorizontalAlignment', 'right');
hs.dpi = uicontrol(d, 'Style', 'popup', 'String', str, ...
    'Position', [110 50 50 22], 'Value', val, 'BackgroundColor', 'w', ...
    'TooltipString', 'in DPI (0 means screen resolution)');
 
uicontrol(d, 'Position', [300 10 70 24], 'Tag', 'OK', ...
    'String', 'OK', 'Callback', @pref_dialog_cb);
uicontrol(d, 'Position',[200 10 70 24], ...
    'String', 'Cancel', 'Callback', 'delete(gcf)');
 
set(d, 'Position', pos, 'Visible', 'on');
guidata(d, hs);
 
%% Preference dialog callback
function pref_dialog_cb(h, ~)
hs = guidata(h);
if h.Tag == "OK" % done
    fh = hs.fig;
    pf = getpref('nii_viewer_para');
    
    pf.layout = hs.layout.Value; % 1:3
    pf.rightOnLeft = hs.rightOnLeft.Value==2;    
    pf.mouseOver = hs.mouseOver.Value;
    if pf.mouseOver % this is the only one we update current fig
        fh.WindowButtonMotionFcn = {@nii_viewer_cb "mousemove" fh};
    else
        fh.WindowButtonMotionFcn = '';        
    end
    
    pf.interp = hs.interp.String{hs.interp.Value};    
    pf.extraV = str2double(hs.extraV.String);
    pf.openPath = hs.openPath.String;
    pf.dpi = hs.dpi.String{hs.dpi.Value};
    
    hs1 = guidata(fh);
    hs1.pref.UserData = pf;
        
    setpref('nii_viewer_para', fieldnames(pf), struct2cell(pf));
    delete(h.Parent);
elseif h.Tag == "browse" % set openPath
    pth = uigetdir(pwd, 'Select folder for background image');
    if ~ischar(pth), return; end
    hs.openPath.String = pth;
end
 
%% Simple version of interp3
function V = interp3a(V, I, method)
% V = interp3a(V, I, 'linear');
% This is similar to interp3 from Matlab, while the input is simplified for
% coordinate. The coordinate input are in this way: I(1,:), I(2,:) and I(3,:)
% are for x, y and z.
if strcmp(method, 'nearest') || any(size(V)<2), I = round(I); end
if size(V,1)<2, V(2,:,:) = nan; end
if size(V,2)<2, V(:,2,:) = nan; end
if size(V,3)<2, V(:,:,2) = nan; end
F = griddedInterpolant(V, method, 'none');
V = F(I(1,:), I(2,:), I(3,:)); % interpolate
 
%% 2D/3D smooth wrapper: no input check for 2D
function im = smooth23(im, varargin)
% out = smooth23(in, varargin)
% This works the same as smooth3 from Matlab, but takes care of 2D input.
is2D = size(im,3) == 1;
if is2D, im = repmat(im, [1 1 2]); end
im = smooth3(im, varargin{:});
if is2D, im = im(:,:,1); end
 
%% Show xyz and value
function xyz = set_xyz(hs, I)
if nargin<2, I = [hs.ijk.Value]; end
I = round(I);
xyz = round(hs.bg.R * [I-1 1]');
xyz = xyz(1:3);
str = sprintf('(%i,%i,%i): ', xyz);
 
for i = 1:size(hs.files.Data,1) % show top one first
    p = get_para(hs, i);
    if p.show == 0, continue; end
    t = round(p.volume);
    if isfield(p, 'R0') % need interpolation
        if isfield(p, 'warp')
            warp = p.warp(I(1), I(2), I(3), :);
            I0 = p.Ri * (p.R0 * [I-1 1]' + [warp(:); 0]);
        else
            I0 = p.Ri * p.R0 * [I-1 1]'; % overlay ijk
        end
        I0 = round(I0(1:3)+1);
    else, I0 = I;
    end
    try
        val = p.nii.img(I0(1), I0(2), I0(3), t, :);
        if isfield(p, 'scl_slope')
            val = single(val) * p.scl_slope + p.scl_inter;
        end
    catch
        val = nan; % out of range
    end
    
    if isfield(p, 'labels')
        try 
            labl = p.labels{val};
            str = sprintf('%s %s', str, labl);
            continue; % skip numeric val assignment
        end
    end
    if isfield(p.nii, 'NamedMap')
        try
            labl = p.nii.NamedMap{p.volume}.labels{val};
            str = sprintf('%s %s', str, labl);
            continue;
        end
    end
    
    fmtstr = '%.4g ';
    if numel(val)>1
        fmtstr = repmat(fmtstr, 1, numel(val));
        fmtstr = ['[' fmtstr]; fmtstr(end) = ']'; %#ok
    end
    str = sprintf(['%s ' fmtstr], str, val);
end
hs.value.Text = str;
 
%% nii essentials
function s = nii_essential(hdr)
% info = nii_essential(hdr);
% Decode important NIfTI hdr into struct with human readable info.
if isfield(hdr, 'nii') % input by nii_viewer
    s.FileName = hdr.nii.hdr.file_name;
    if isfield(hdr, 'mask_info'), s.MaskedBy = hdr.mask_info; end
    if isfield(hdr, 'modulation_info'), s.ModulatedBy = hdr.modulation_info; end
    if isfield(hdr, 'alignMtx'), s.AlignMatrix = hdr.alignMtx; end
    hdr = hdr.hdr0; % original nii header
else
    s.FileName = hdr.file_name;
end
switch hdr.intent_code
    case 2, s.intent = 'Correlation'; s.DoF = hdr.intent_p1;
    case 3, s.intent = 'T-test';      s.DoF = hdr.intent_p1;
    case 4, s.intent = 'F-test';      s.DoF = [hdr.intent_p1 hdr.intent_p2];
    case 5, s.intent = 'Z-score';
    case 6, s.intent = 'Chi squared'; s.DoF = hdr.intent_p1;
        
    % several non-statistical intent_code
    case 1002, s.intent = 'Label'; % e.g. AAL labels
    case 2003, s.intent = 'RGB'; % triplet in the 5th dimension
    case 2004, s.intent = 'RGBA'; % quadruplet in the 5th dimension
end
switch hdr.datatype
    case 0
    case 1,    s.DataType = 'logical';
    case 2,    s.DataType = 'uint8';
    case 4,    s.DataType = 'int16';
    case 8,    s.DataType = 'int32';
    case 16,   s.DataType = 'single';
    case 32,   s.DataType = 'single complex';
    case 64,   s.DataType = 'double';
    case 128,  s.DataType = 'uint8 RGB';
    case 256,  s.DataType = 'int8';
    case 511,  s.DataType = 'single RGB';
    case 512,  s.DataType = 'uint16';
    case 768,  s.DataType = 'uint32';
    case 1024, s.DataType = 'int64';
    case 1280, s.DataType = 'uint64';
    case 1792, s.DataType = 'double complex';
    case 2304, s.DataType = 'uint8 RGBA';
    otherwise, s.DataType = 'unknown';
end
s.Dimension = hdr.dim(2:hdr.dim(1)+1);
switch bitand(hdr.xyzt_units, 7)
    case 1, s.VoxelUnit = 'meters';
    case 2, s.VoxelUnit = 'millimeters';
    case 3, s.VoxelUnit = 'micrometers';
end 
s.VoxelSize = hdr.pixdim(2:4);
switch bitand(hdr.xyzt_units, 56)
    case 8,  s.TemporalUnit = 'seconds';
    case 16, s.TemporalUnit = 'milliseconds';
    case 24, s.TemporalUnit = 'microseconds';
    case 32, s.TemporalUnit = 'Hertz';
    case 40, s.TemporalUnit = 'ppm';
    case 48, s.TemporalUnit = 'radians per second';
end 
if isfield(s, 'TemporalUnit') && strfind(s.TemporalUnit, 'seconds')
    s.TR = hdr.pixdim(5);
end
if hdr.dim_info>0
    s.FreqPhaseSliceDim = bitand(hdr.dim_info, [3 12 48]) ./ [1 4 16];
    a = bitand(hdr.dim_info, 192) / 64;
    if a>0 && s.FreqPhaseSliceDim(2)>0
        ax = 'xyz'; % ijk to be accurate
        pm = ''; if a == 2, pm = '-'; end 
        s.PhaseEncodingDirection = [pm ax(s.FreqPhaseSliceDim(2))]; 
    end
end
 
switch hdr.slice_code
    case 0
    case 1, s.SliceOrder = 'Sequential Increasing';
    case 2, s.SliceOrder = 'Sequential Decreasing';
    case 3, s.SliceOrder = 'Alternating Increasing 1';
    case 4, s.SliceOrder = 'Alternating Decreasing 1';
    case 5, s.SliceOrder = 'Alternating Increasing 2';
    case 6, s.SliceOrder = 'Alternating Decreasing 2';
    otherwise, s.SliceOrder = 'Multiband?';
end
if ~isempty(hdr.descrip), s.Notes = hdr.descrip; end
str = formcode2str(hdr.qform_code);
if ~isempty(str), s.qform = str; end
if hdr.qform_code>0, s.qform_mat = hdr.qform_mat; end
str = formcode2str(hdr.sform_code);
if ~isempty(str), s.sform = str; end
if hdr.sform_code>0, s.sform_mat = hdr.sform_mat; end
 
%% decode NIfTI form_code
function str = formcode2str(code)
switch code
    case 0, str = 'Undefined';
    case 1, str = 'Scanner Anat';
    case 2, str = 'Aligned Anat';
    case 3, str = 'Talairach';
    case 4, str = 'mni_152';
    otherwise, str = sprintf('Formcode %i', code);
end
 
%% Get a mask based on image intensity, but with inside brain filled
function r = img2mask(img, thr)
if nargin<2 || isempty(thr), thr = mean(img(img(:)>0)) / 8; end
r = smooth23(img, 'box', 5) > thr; % smooth, binarize
if sum(r(:))==0, return; end
C = contourc(double(r), [1 1]);
i = 1; c = {};
while size(C,2)>2 % split C into contours
    k = C(2,1) + 1;
    c{i} = C(:, 2:k); C(:,1:k) = []; %#ok
    i = i+1;
end

nc = numel(c);
rg = nan(nc, 4); % minX minY maxX maxY
for i = 1:nc
    rg(i,:) = [min(c{i},[],2)' max(c{i},[],2)'];
end
ind = false(nc,1);
foo = min(rg(:,1)); ind = ind | foo==rg(:,1);
foo = min(rg(:,2)); ind = ind | foo==rg(:,2);
foo = max(rg(:,3)); ind = ind | foo==rg(:,3);
foo = max(rg(:,4)); ind = ind | foo==rg(:,4);
c = c(ind); % outmost contour(s)
len = cellfun(@(x) size(x,2), c);
[~, ind] = sort(len, 'descend');
c = c(ind);
C = c{1};
if isequal(C(:,1), C(:,end)), c(2:end) = []; end % only 1st if closed
nc = numel(c);
for i = nc:-1:2 % remove closed contours except one with max len
    if isequal(c{i}(:,1), c{i}(:,end)), c(i) = []; end
end
nc = numel(c);
while nc>1 % +1 contours, put all into 1st
    d2 = nan(nc-1, 2); % distance^2 from C(:,end) to other start/endpoint
    for i = 2:nc
        d2(i-1,:) = sum((C(:,end)*[1 1] - c{i}(:,[1 end])).^2);
    end
    [i, j] = find(d2 == min(d2(:)));
    i = i + 1; % start with 2nd
    if j == 1, C = [C c{i}]; %#ok C(:,1) connect to c{i}(:,1)
    else C = [C c{i}(:,end:-1:1)]; %#ok C(:,end) to c{i}(:,end)
    end
    c(i) = []; nc = nc-1;
end
if ~isequal(C(:,1), C(:,end)), C(:,end+1) = C(:,1); end % close the contour
x = C(1, :);
y = C(2, :);
[m, n] = size(r);

% following is the method in Octave poly2mask
xe = [x(1:numel(x)-1); x(1, 2:numel(x))]; % edge x
ye = [y(1:numel(y)-1); y(1, 2:numel(y))]; % edge y
ind = ye(1,:) == ye(2,:);
xe(:,ind) = []; ye(:, ind) = []; % reomve horizontal edges
minye = min(ye);
maxye = max(ye);
t = (ye == [minye; minye]);
exminy = xe(:); exminy = exminy(t);
exmaxy = xe(:); exmaxy = exmaxy(~t);
maxye = maxye';
minye = minye';
m_inv = (exmaxy - exminy) ./ (maxye - minye);
ge = [maxye minye exmaxy m_inv];
ge = sortrows(ge, [1 3]);
ge = [-Inf -Inf -Inf -Inf; ge];

gei = size(ge, 1);
sl = ge(gei, 1);
ae = []; % active edge
while (sl == ge(gei, 1))
    ae = [ge(gei, 2:4); ae]; %#ok
    gei = gei - 1;
end

miny = min(y);
if miny < 1, miny = 1; end

while (sl >= miny)
    if (sl <= m) % check vert clipping
        ie = round(reshape(ae(:, 2), 2, size(ae, 1)/2));
        ie(1, :) = ie(1, :) + (ie(1, :) ~= ie(2, :));
        ie(1, (ie(1, :) < 1)) = 1;
        ie(2, (ie(2, :) > n)) = n;
        ie = ie(:, (ie(1, :) <= n));
        ie = ie(:, (ie(2, :) >= 1));
        for i = 1:size(ie,2)
            r(sl, ie(1, i):ie(2, i)) = true;
        end
    end

    sl = sl - 1;
    ae = ae((ae(:, 1) ~= sl), :);
    ae(:, 2) = ae(:, 2) - ae(:, 3);

    while (sl == ge(gei, 1))
        ae = [ae; ge(gei, 2:4)]; %#ok
        gei = gei - 1;
    end

    if size(ae,1) > 0
        ae = sortrows(ae, 2);
    end
end
 
%% update colorbar label
function set_colorbar(hs)
if strcmpi(hs.colorbar.Visible, 'off'), return; end
p = get_para(hs);
map = lut2map(p);
rg = sort([p.lb p.ub]);
tickLoc = [0 0.5 1];
if startsWith(p.lut, "phase")
    labls = [0 180 360];
elseif p.lut == "two-sided"
    rg = sort(abs(rg));
    labls = {num2str(-rg(2),'%.3g') num2str(rg(1),'+/-%.3g') num2str(rg(2),'%.3g')};
elseif p.lut == "custom"
    im = p.nii.img(:,:,:,p.volume);
    im(isnan(im) | im==0) = [];
    im = unique(im);
    if max(im)<=size(map,1) && isequal(im, round(im)) % integer
        try map = map(im+1, :); rg = [im(1) im(end)]; end
    end
    labls = {num2str(rg(1),'%.2g') num2str(rg(2),'%.3g')};
    tickLoc = [0 1];
else
    if rg(2)<0, rg = rg([2 1]); end
    mn = str2double(num2str(mean(rg), '%.4g'));
    labls = [rg(1) mn rg(2)];
end
colormap(hs.ax(end), map); % map must be double for old matlab
set(hs.colorbar, 'YTickLabel', labls, 'YTick', tickLoc, 'Ylim', [0 1]);
 
fh = hs.fig.UserData;
if isempty(fh) || ~ishandle(fh) || ~isfield(p.nii, 'cii'), return; end
hs = guidata(fh);
colormap(hs.ax(end), map);
set(hs.colorbar, 'YTickLabel', labls, 'YTick', tickLoc, 'Ylim', [0 1]);
 
%% return screen size in pixels
function res = screen_pixels(id)
res = get(0, 'MonitorPositions');
if size(res,1)<2, res = res(1, 3:4); return; end % single/duplicate monitor
if nargin, res = res(id,3:4); return; end
res = sortrows(res);
res = res(end,1:2) + res(end,3:4) - res(1,1:2);
 
%% add mask or modulation
function addMask(h, ~)
hs = guidata(h);
jf = -hs.stack.Value;
p = get_para(hs, jf);
pName = fileparts(p.nii.hdr.file_name);
[fname, pName] = uigetfile([pName '/*.nii;*.hdr;*.nii.gz;*.hdr.gz'], 'Select mask NIfTI');
if ~ischar(fname), return; end
fname = fullfile(pName, fname);
 
nii = nii_tool('load', fname);
hdr = p.nii.hdr;
codes = [hdr.sform_code hdr.qform_code];
[R, frm] = nii_xform_mat(nii.hdr, codes);
if ~any(frm == codes)
    str = ['There is no matching coordinate system between the selected ' ...
        'image and the mask image. Do you want to apply the mask anyway?'];
    btn = questdlg(str, 'Apply mask?', 'Cancel', 'Apply', 'Cancel');
    if isempty(btn) || strcmp(btn, 'Cancel'), return; end
    R0 = nii_xform_mat(hdr, codes(1));
else
    R0 = nii_xform_mat(hdr, frm); % may be the same as p.R
end
% R0 = reorient(R0, hdr.dim(2:4)); % do this since img was done when loaded
 
% if isfield(p, 'alignMtx'), R = R0 / p.R * R; end % inverse
% this wont work if lines is background & Mprage is overlay
if all(isfield(p, {'R0' 'alignMtx'})) % target as mask
    R1 = reorient(R, nii.hdr.dim(2:4));
    if all(abs(R1(:)-p.R0(:))<1e-4), R0 = p.R; end % not 100% safe
end
 
d = single(size(p.nii.img)); % dim for reoriented img
d(numel(d)+1:3) = 1; d = d(1:3);
 
I = inv(R) * R0 * grid3(d) + 1; %#ok ijk+1 for mask
I = round(I * 100) / 100;
 
im = single(nii.img(:,:,:,1)); % first mask volume
slope = nii.hdr.scl_slope;
if slope==0, slope = 1; end
im = im * slope + nii.hdr.scl_inter;
im = interp3a(im, I, 'nearest');
im1 = im(~isnan(im)); % for threshold computation
im = reshape(im, d);
 
if h.Label == "Apply mask" % binary mask
    if numel(unique(im1))<3
        thre = min(im1);
    else
        a = get_range(nii);
        str = sprintf('Threshold for non-binary mask (%.3g to %.4g)', ...
            min(im1), max(im1));
        a = inputdlg(str, 'Input mask threshold', 1, {num2str(a(1), '%.3g')});
        if isempty(a), return; end
        thre = str2double(a{1});
        fname = [fname ' (threshold = ' a{1} ')']; % info only
    end
    p.mask = ones(size(im), 'single');
    p.mask(abs(im)<=thre) = nan;
    p.mask_info = fname;
    noteStr = '(masked)';
else % modulation
    mi = min(im1); ma = max(im1);
    if mi<0 || ma>1
        str = {sprintf('Lower bound to clip to 0 (image min = %.2g)', mi) ...
               sprintf('Upper bound to clip to 1 (image max = %.2g)', ma)};
        def = strtrim(cellstr(num2str([mi;ma], '%.2g'))');
        a = inputdlg(str, 'Input modulation range', 1, def);
        if isempty(a), return; end
        mi = str2double(a{1}); ma = str2double(a{2});
        fname = [fname ' (range [' a{1} ' ' a{2} '])'];
    end
    im(im<mi) = mi; im(im>ma) = ma;
    p.modulation = (im-mi) / (ma-mi);
    p.modulation_info = fname;
    noteStr = '(modulated)';
end
p.hsI(1).UserData = p;
set_cdata(hs);
 
str = hs.files.Data{jf,2};
if ~any(regexp(str, [regexptranslate('escape', noteStr) '$']))
    hs.files.Data{jf,2} = [str noteStr];
    table_width(hs.files);
end
 
%% Return 0-based 4xN 3D grid: [i; j; k; 1]
function I = grid3(d)
I = ones([4 d(1:3)], 'single');
[I(1,:,:,:), I(2,:,:,:), I(3,:,:,:)] = ndgrid(0:d(1)-1, 0:d(2)-1, 0:d(3)-1);
I = reshape(I, 4, []);
 
%% update crosshair: ix correspond to one of the three spinners, not views
function set_cross(hs, ix)
h = hs.cross;
for i = ix
    c = hs.ijk(i).Value;
    g = hs.gap(i);
    if i == 1 % I changed
        set([h(2,3:4) h(3,3:4)], 'XData', [c c]);
        set([h(2,1) h(3,1)], 'XData', [0 c-g]);
        set([h(2,2) h(3,2)], 'XData', [c+g hs.dim(1)+1]);
    elseif i == 2 % J
        set(h(1,3:4), 'XData', [c c]);
        set(h(1,1), 'XData', [0 c-g]);
        set(h(1,2), 'XData', [c+g hs.dim(2)+1]);
        set(h(3,1:2), 'YData', [c c]);
        set(h(3,3), 'YData', [0 c-g]);
        set(h(3,4), 'YData', [c+g hs.dim(2)+1]);
    else % K
        set([h(1,1:2) h(2,1:2)], 'YData', [c c]);
        set([h(1,3) h(2,3)], 'YData', [0 c-g]);
        set([h(1,4) h(2,4)], 'YData', [c+g hs.dim(3)+1]);
    end
end
 
%% Duplicate image handles, including ButtonDownFcn for new matlab
function h = copyimg(hs)
h = hs.hsI;
cb = {@nii_viewer_cb 'mousedown' hs.fig};
if isvalid(hs.hsI(1))
    for i = 1:3, h(i) = handle(copyobj(hs.hsI(i), hs.ax(i))); end
    set(h, 'Visible', 'on', 'ButtonDownFcn', cb);
else % when lut=lines is the background
    for i = 1:3
        j = 1:3; j(j==i) = [];
        hs.hsI(i) = handle(image(hs.ax(i), zeros(hs.dim(j([2 1])), 'single')));
    end
    set(hs.hsI, 'Visible', 'on', 'ButtonDownFcn', cb);
    guidata(hs.fig, hs);
    h = hs.hsI;
end
crossFront(hs);
 
%% Save selected nii as another file
function save_nii_as(h, ~)
fh = ancestor(h, 'figure');
hs = guidata(fh);
labl = h.Label;
p = get_para(hs);
nam = p.nii.hdr.file_name;
pName = fileparts(nam);
if isempty(pName), pName = pwd; end
try nii = nii_tool('load', nam); % re-load to be safe
catch % restore reoriented img
    nii = p.nii; % just save reoriented version if no original img
    slope = nii.hdr.scl_slope; if slope==0, slope = 1; end
    nii.img = (single(nii.img) - nii.hdr.scl_inter) / slope; % undo scale
    if nii.hdr.datatype == 4 % leave others as it is or single
        nii.img = int16(nii.img);
    elseif nii.hdr.datatype == 512
        nii.img = uint16(nii.img);
    elseif any(nii.hdr.datatype == [2 128 2304])
        nii.img = uint8(nii.img);
    end
end
 
if contains(labl, 'a copy') % a copy
    [fname, pName] = uiputfile([pName '/*.nii'], 'Input file name');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    nii_tool('save', nii, fname);
elseif contains(labl, 'dim 4') % fsl RGB
    if any(size(nii.img,8) == 3:4)
        nii.img = permute(nii.img, [1:3 8 4:7]);
    elseif ~any(nii.hdr.dim(5) == 3:4)
        uialert(fh, 'Selected image is not RGB data.', 'Error'); return;
    end
    [fname, pName] = uiputfile([pName '/*.nii'], 'Input name for FSL RGB file');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    nii_tool('save', nii, fname);
elseif contains(labl, 'dim 3') % old mricron RGB
    if any(nii.hdr.dim(5) == 3:4)
        nii.img = permute(nii.img, [1:3 5:7 4]);
    elseif ~any(size(nii.img,8) == 3:4)
        uialert(fh, 'Selected image is not RGB data', 'Error'); return;
    end
    [fname, pName] = uiputfile([pName '/*.nii'], 'Input name for old mricrom styte file');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    old = nii_tool('RGBStyle', 'mricron');
    nii_tool('save', nii, fname);
    nii_tool('RGBStyle', old);
elseif contains(labl, 'AFNI') % NIfTI RGB
    if any(nii.hdr.dim(5) == 3:4)
        nii.img = permute(nii.img, [1:3 5:8 4]);
    elseif ~any(size(nii.img,8) == 3:4)
        uialert(fh, 'Selected image is not RGB data', 'Error'); return;
    end
    nii.img = abs(nii.img);
    [fname, pName] = uiputfile([pName '/*.nii'], ...
        'Input name for NIfTI standard RGB file');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    old = nii_tool('RGBStyle', 'afni');
    nii_tool('save', nii, fname);
    nii_tool('RGBStyle', old);
elseif contains(labl, '3D') % SPM 3D
    if nii.hdr.dim(5)<2
        uialert(fh, 'Selected image is not multi-volume data', 'Error'); return;
    end
    [fname, pName] = uiputfile([pName '/*.nii'], 'Input base name for SPM 3D file');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    nii_tool('save', nii, fname, 1); % force 3D
elseif contains(labl, 'new resolution')
    str = 'Resolution for three dimension in mm';
    a = inputdlg(str, 'Input spatial resolution', 1, {'3 3 3'});
    if isempty(a), return; end
    res = sscanf(a{1}, '%g %g %g');
    if numel(res) ~= 3
        uialert(fh, 'Invalid spatial resolution', 'Error');
        return;
    end
    if isequal(res, nii.hdr.pixdim(2:4))
        uialert(fh, 'The input resolution is the same as current one', 'Error');
        return;
    end
    [fname, pName] = uiputfile([pName '/*.nii;nii.gz'], ...
        'Input result name for the new resolution file');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    nii_xform(nii, res, fname, hs.pref.UserData.interp, hs.pref.UserData.extraV)
elseif contains(labl, 'matching background')
    if p.hsI(1) == hs.hsI(1)
        uialert(fh, 'You selected background image', 'Error');
        return;
    end
    [fname, pName] = uiputfile([pName '/*.nii;*.nii.gz'], ...
        'Input result file name');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    nii_xform(nii, hs.bg.hdr, fname, hs.pref.UserData.interp, hs.pref.UserData.extraV)
elseif contains(labl, 'aligned template')
    [temp, pName] = uigetfile([pName '/*.nii;*.nii.gz'], ...
        'Select the aligned template file');
    if ~ischar(temp), return; end
    temp = fullfile(pName, temp);
    [mtx, pName] = uigetfile([pName '/*.mat'], ['Select the text ' ...
        'matrix file which aligns the nii to the template']);
    if ~ischar(mtx), return; end
    mtx = fullfile(pName, mtx);
    [fname, pName] = uiputfile([pName '/*.nii;*.nii.gz'], ...
        'Input result file name');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    nii_xform(nii, {temp mtx}, fname, hs.pref.UserData.interp, hs.pref.UserData.extraV)
else
    uialert(fh, sprintf('%s not implemented yet.', labl), 'Warning');
end
 
%% Return 3-layer RGB, called by set_cdata
function [im, alfa] = lut2img(im, p)
rg = sort([p.lb p.ub]);
try p.map = p.nii.NamedMap{p.volume}.map; end
if startsWith(p.lut, "phase")
    im = im(:,:,2) .* single(im(:,:,1)>min(abs(rg))); % mag as mask
end
if rg(2)<0 % asking for negative data
    rg = -rg([2 1]);
    if p.lut~="two-sided", im = -im; end
end
 
alfa = single(0); % override for lut=10
if p.lut == "two-sided" % store negative value
    rg = sort(abs(rg));
    im_neg = -single(im) .* (im<0);
    im_neg = (im_neg-rg(1)) / (rg(2)-rg(1));
    im_neg(im_neg>1) = 1; im_neg(im_neg<0) = 0;
    alfa = im_neg; % add positive part later
    im_neg = repmat(im_neg, [1 1 3]); % gray now
end
 
if ~any([startsWith(p.lut, "phase") (p.lut == ["RGB" "custom"])]) % no scaling
    im = (im-rg(1)) / (rg(2)-rg(1));
    im(im>1) = 1; im(im<0) = 0;
end
alfa = im + alfa;
if p.lut ~= "RGB", im = repmat(im, [1 1 3]); end % gray now
 
switch p.lut
    case "grayscale" % do nothing
    case "red",     im(:,:,2:3) = 0;
    case "green",   im(:,:,[1 3]) = 0;
    case "blue",    im(:,:,1:2) = 0;
    case "violet",  im(:,:,2) = 0;
    case "yellow",  im(:,:,3) = 0;
    case "cyan",    im(:,:,1) = 0;
    case {"red-yellow" "autumn"}
        a = im(:,:,1); a(a>0) = 1; im(:,:,1) = a;
        im(:,:,3) = 0;
    case "blue-green"
        im(:,:,1) = 0;
        a = im(:,:,3); a(a==0) = 1; a = 1 - a; im(:,:,3) = a;
    case "two-sided" % combine red_yellow & blue_green
        im(:,:,3) = 0;
        a = im(:,:,1); a(a>0) = 1; im(:,:,1) = a;
        im_neg(:,:,1) = 0;
        a = im_neg(:,:,3); a(a==0) = 1; a = 1 - a; im_neg(:,:,3) = a;
        im = im + im_neg;
    case "hot" % Matlab colormap can be omitted, but faster than mapping
        a = im(:,:,1); a = a/0.375; a(a>1) = 1; im(:,:,1) = a;
        a = im(:,:,2); a = a/0.375-1;
        a(a<0) = 0; a(a>1) = 1; im(:,:,2) = a;
        a = im(:,:,3); a = a*4-3; a(a<0) = 0; im(:,:,3) = a;
    case "cool"
        a = im(:,:,2); a(a==0) = 1; a = 1 - a; im(:,:,2) = a;
        a = im(:,:,3); a(a>0) = 1; im(:,:,3) = a;
    case "spring"
        a = im(:,:,1); a(a>0) = 1; im(:,:,1) = a;
        a = im(:,:,3); a(a==0) = 1; a = 1 - a; im(:,:,3) = a;
    case "summer"
        a = im(:,:,2); a(a==0) = -1; a = a/2+0.5; im(:,:,2) = a;
        a = im(:,:,3); a(a>0) = 0.4; im(:,:,3) = a;
    case "winter"
        im(:,:,1) = 0;
        a = im(:,:,3); a(a==0) = 2; a = 1-a/2; im(:,:,3) = a;
    case "copper"
        a = im(:,:,1); a = a*1.25; a(a>1) = 1; im(:,:,1) = a;
        im(:,:,2) = im(:,:,2) * 0.7812;
        im(:,:,3) = im(:,:,3) * 0.5;
    case "phase" % like red_yellow
        im(:,:,1) = 1; im(:,:,3) = 0;
    case "phase3" % red-yellow-green-yellow-red
        a = im(:,:,1);
        b1 = a<=0.25;
        b2 = a>0.25 & a<=0.5;
        b3 = a>0.5 & a<=0.75;
        b4 = a>0.75;
        a(b1 | b4) = 1;
        a(b2) = (0.5-a(b2))*4;
        a(b3) = (a(b3)-0.5)*4;
        im(:,:,1) = a;
        
        a = im(:,:,2);
        a(b2 | b3) = 1;
        a(b1) = a(b1)*4;
        a(b4) = (1-a(b4))*4;
        im(:,:,2) = a;
        
        im(:,:,3) = 0;
    case "phase6" % red-yellow-green/violet-blue-cyan
        a = im(:,:,1);
        b1 = a<=0.25;
        b2 = a>0.25 & a<=0.5;
        b3 = a>0.5 & a<=0.75;
        b4 = a>0.75;
        a(b2) =  (0.5-a(b2))*4; a(b1) = 1;
        a(b3) = (0.75-a(b3))*4; a(b4) = 0;
        im(:,:,1) = a;
        
        a = im(:,:,2);
        a(b1) = a(b1)*4; a(b2) = 1;
        a(b3) = 0; a(b4) = (a(b4)-0.75)*4;
        im(:,:,2) = a;
        
        a = im(:,:,3);
        a(b1 | b2) = 0;
        a(b3 | b4) = 1;
        im(:,:,3) = a;
    case "RGB" % disp non-NIfTI RGB as RGB
        im = abs(im); % e.g. DTI V1
        if max(im(:)) > 1, im = im / 255; end % it should be unit8
        alfa = sum(im,3)/3;
    otherwise % parula(12), jet(13), hsv(14), bone(21), pink(23), custom
        if isfield(p, 'map') % custom
            map = p.map;            
        else
            map = feval(p.lut, 256);
        end
        if p.lut ~= "custom" % normalized previously
            a = floor(im(:,:,1) * (size(map,1)-1)) + 1; % 1st for bkgrnd
        elseif max(p.nii.img(:)) <= size(map,1)
            alfa = alfa / max(alfa(:));
            a = round(im(:,:,1)) + 1; % custom or uint8, round to be safe
        else
            a = (im(:,:,1) - rg(1)) / (rg(2)-rg(1));
            a(a<0) = 0 ; a(a>1) = 1;
            alfa = a;
            a = round(a * (size(map,1)-1)) + 1;
        end
        a(isnan(a)) = 1;
        aa = a;
        for i = 1:3
            aa(:) = map(a, i);
            im(:,:,i) = aa;
        end
end
 
%% Return binary sphere ROI from xyz and r (mm)
function b = xyzr2roi(c, r, hdr)
% ROI_img = xyzr2roi(center, radius, hdr)
% Return an ROI img based on the dim info in NIfTI hdr. The center and radius
% are in unit of mm. 
d = single(hdr.dim(2:4));
I = nii_xform_mat(hdr) * grid3(d); % xyz in 4 by nVox
 
b = I(1:3,:) - c(:); % dist in x y z direction from center
b = sum(b .* b); % dist to center squared, 1 by nVox
 
b = b <= r*r; % within sphere
b = reshape(b, d);
 
%% Return center of gravity of an image
function c = img_cog(img)
% center_ijk = img_cog(img)
% Return the index of center of gravity in img (must be 3D).
img(isnan(img)) = 0;
img = double(abs(img));
img = img / sum(img(:));
c = ones(3,1);
for i = 1:3
    a = permute(img, [i 1:i-1 i+1:3]);
    c(i) = (1:size(img,i)) * sum(sum(a,3),2);
end
 
%% set up disp parameter for new nifti
function p = dispPara(p, luts)
if nargin<2, luts = []; end
p.show = true; % img on
if isfield(p, 'map')
    p.lut = "custom";
elseif any(p.nii.hdr.datatype == [32 1792]) % complex
    p.lut = "phase";
    p.lb = str2double(sprintf('%.2g', p.ub/2));
elseif any(p.nii.hdr.intent_code == [1002 3007]) % Label
    p.lut = "prism";
elseif size(p.nii.img,8)>1
    p.lut = "RGB";
elseif p.nii.hdr.intent_code > 0 % some stats
    if p.lb < 0
        p.lut = "two-sided";
        p.lb = round(p.ub/2, 2, 'significant');
    else
        a = setdiff(["red-yellow" "blue-green"], luts, 'stable');
        if isempty(a), a = "red-yellow"; end
        p.lut = a(1);
    end
elseif isempty(luts)
    p.lut = "grayscale";
else
    a = setdiff(["red" "green" "blue" "violet" "yellow" "cyan"], luts, 'stable');
    if isempty(a), a = "red"; end
    p.lut = a(1);
end
p.lb_step = stepSize(p.lb); 
p.ub_step = stepSize(p.ub);
p.alpha = 1; % opaque
p.smooth = false;
p.interp = "nearest";
p.volume = 1; % first volume
 
%% estimate StepSize for spinner
function d = stepSize(val)
d = round(abs(val/10), 1, 'significant');
d = max(d, 0.01);
if d>4, d = round(d/2)*2; end
 
%% Return nii struct from nii struct, nii fname or other convertible files
function nii = get_nii(fname)
if isstruct(fname), nii = fname; return;
elseif iscellstr(fname), nam = fname{1}; %#ok
else, nam = fname;
end
try 
    nii = nii_tool('load', strtrim(nam));
catch me
    try nii = dicm2nii(fname, pwd, 'no_save');
    catch, rethrow(me);
    end
end
 
%% Get figure/plot position from FoV for layout
% siz is in pixels, while pos is normalized.
function [siz, axPos, figPos] = plot_pos(mm, layout)
if layout==1 % 1x3
    siz = [sum(mm([2 1 1]))+mm(1)/4 max(mm(2:3))]; % image area width/height
    y0 = mm(2) / siz(1); % normalized width of sag images
    x0 = mm(1) / siz(1); % normalized width of cor/tra image
    z0 = mm(3) / siz(2); % normalized height of sag/cor images
    y1 = mm(2) / siz(2); % normalized height of tra image
    if y1>z0, y3 = 0; y12 = (y1-z0)/2;
    else, y3 = (z0-y1)/2; y12 = 0;
    end
    axPos = [0 y12 y0 z0;  y0 y12 x0 z0;  y0+x0 y3 x0 y1;  y0+x0*2 0 mm(1)/4/siz(1) min(z0,y1)];
elseif layout==2 || layout==3 % 2x2
    siz = [sum(mm(1:2)) sum(mm(2:3))]; % image area width/height
    x0 = mm(1) / siz(1); % normalized width of cor/tra images
    y0 = mm(2) / siz(2); % normalized height of tra image
    z0 = mm(3) / siz(2); % normalized height of sag/cor images
    if layout == 2 % 2x2 sag at (1,2)
        axPos = [x0 y0 1-x0 z0;  0 y0 x0 z0;  0 0 x0 y0;  x0 0 1-x0 y0];
    else % ==3:      2x2 sag at (1,2)
        axPos = [0 y0 1-x0 z0;  1-x0 y0 x0 z0;  1-x0 0 x0 y0;  0 0 1-x0 y0];
    end
else
    error('Unknown layout parameter');
end
siz = siz / max(siz) * 800;
 
res = screen_pixels(1); % use 1st screen
maxH = res(2) - 180;
maxW = res(1) - 100;
if siz(1)>maxW, siz = siz / siz(1) * maxW; end
if siz(2)>maxH, siz = siz / siz(2) * maxH; end
 
figPos = round((res-siz)/2);
if figPos(1)+siz(1) > res(1), figPos(1) = res(1)-siz(1)-10; end
if figPos(2)+siz(2) > res(2)-180, figPos(2) = min(figPos(2), 50); end
 
%% Return nii struct from cii and gii
function nii = cii2nii(nii)
persistent gii; % Anatomical surface
if nargin<1, nii = gii; return; end
 
ind = find([nii.ext.ecode]==32, 1);
xml = nii.ext(ind).edata_decoded;
expr = '(?<=((SurfaceNumberOfVertices)|(SurfaceNumberOfNodes))=").*?(?=")';
nVer = str2double(regexp(xml, expr, 'match', 'once'));
if isempty(nVer), error('SurfaceNumberOfVertices not found'); end
if isempty(gii), gii = get_surfaces(nVer, 'Anatomical'); end
if isempty(gii) || numel(gii.Vertices) ~=2, error('Not valid GIfTI'); end
if nVer ~= size(gii.Vertices{1},1), error('GIfTI and CIfTI don''t match'); end
 
if gii_attr(xml, 'CIFTI Version', 1) == 1
    dim = size(nii.img);
    nii.img = reshape(nii.img, dim([1:4 6 5]));
end
 
dim = gii_attr(xml, 'VolumeDimensions', 1);
if isempty(dim), dim = [91 109 91]; end % HCP 2x2x2 mm
TR = gii_attr(xml, 'SeriesStep', 1);
if ~isempty(TR), nii.hdr.pixdim(5) = TR; end
mat = gii_element(xml, 'TransformationMatrixVoxelIndicesIJKtoXYZ', 1);
 
if isempty(mat) % some cii miss 'mat' and 'dim'
    mat = [-2 0 0 90; 0 2 0 -126; 0 0 2 -72; 0 0 0 1]; % best guess from HCP
else
    pow = gii_attr(xml, 'MeterExponent', 1);
    mat(1:3,:) = mat(1:3,:) / 10^(3+pow);
end
nii.hdr.sform_code = gii.DataSpace;
nii.hdr.sform_mat = mat(1:3,:);
nii.hdr.pixdim(2:4) = sqrt(sum(mat(1:3,1:3).^2));
 
nVol = size(nii.img, 5);
imgG = permute(nii.img, [6 5 1:4]);
nii.img = zeros([prod(dim) nVol], class(imgG));
 
iMdl = regexp(xml, '<BrainModel[\s>]'); iMdl(end+1) = numel(xml);
for j = 1:numel(iMdl)-1
    c = xml(iMdl(j):iMdl(j+1));
    offset = gii_attr(c, 'IndexOffset', 1);
    typ = gii_attr(c, 'ModelType');
    if strcmp(typ, 'CIFTI_MODEL_TYPE_SURFACE')
        a = gii_attr(c, 'BrainStructure');
        ig = find(strcmp(a, {'CIFTI_STRUCTURE_CORTEX_LEFT' 'CIFTI_STRUCTURE_CORTEX_RIGHT'}));
        if isempty(ig), warning('Unknown BrainStructure: %s', a); continue; end
        ind = gii_element(c, 'VertexIndices', 1) + 1;
        if isempty(ind), ind = 1:gii_attr(c, 'IndexCount', 1); end
        nii.cii{ig} = zeros(nVer, nVol, 'single');
        nii.cii{ig}(ind,:) = imgG((1:numel(ind))+offset, :);
    elseif strcmp(typ, 'CIFTI_MODEL_TYPE_VOXELS')
        a = gii_element(c, 'VoxelIndicesIJK', 1) + 1;
        a = sub2ind(dim, a(:,1), a(:,2), a(:,3));
        nii.img(a, :) = imgG((1:numel(a))+offset, :);
    end
end
 
for ig = 1:2 % map surface back to volume
    v = gii.Vertices{ig}'; v(4,:) = 1;
    v = round(mat \ v) + 1; % ijk 1-based
    ind = sub2ind(dim, v(1,:), v(2,:), v(3,:));
    nii.img(ind, :) = nii.cii{ig};
end
nii.img = reshape(nii.img, [dim nVol]);
nii = nii_tool('update', nii);
 
expr = '<Label\s+Key="(.*?)"\s+Red="(.*?)"\s+Green="(.*?)"\s+Blue="(.*?)".*?>(.*?)</Label>';
ind = regexp(xml, '<NamedMap'); ind(end+1) = numel(xml);
for k = 1:numel(ind)-1
    nii.NamedMap{k}.MapName = gii_element(xml(ind(k):ind(k+1)), 'MapName');
    tok = regexp(xml(ind(k):ind(k+1)), expr, 'tokens');
    if numel(tok)<2, continue; end
    tok = reshape([tok{:}], 5, [])'; % Key R G B nam
    a = str2double(tok(:, 1:4));
    nii.NamedMap{k}.map(a(:,1)+1, :) = a(:,2:4);
    if a(1,1) == 0 % key="0"
        nii.NamedMap{k}.labels(a(2:end, 1)) = tok(2:end, 5); % drop Key="0"
    else
        nii.NamedMap{k}.map = [0 0 0; nii.NamedMap{k}.map];
        nii.NamedMap{k}.labels(a(:, 1)) = tok(:, 5);
    end
end
 
%% Return gii for both hemesperes
function gii = get_surfaces(nVer, surfType)
persistent pth;
if nargin>1 && strcmpi(surfType, 'Anatomical') && nVer == 32492 % HCP 2mm surface
    fname = fullfile(fileparts(mfilename('fullpath')), 'example_data.mat');
    a = load(fname, 'gii'); gii = a.gii;
    return;
end
if nargin>1, surfType = [surfType ' ']; else, surfType = ''; end
prompt = ['Select ' surfType 'GIfTI surface files for both hemispheres'];
if isempty(pth), pth = pwd; end
[nam, pth] = uigetfile([pth '/*.surf.gii'],  prompt, 'MultiSelect', 'on');
if isnumeric(nam), gii = []; return; end
nam = cellstr(strcat([pth '/'], nam));
for i = 1:numel(nam)
    a = read_gii(nam{i});
    if ~isempty(surfType) && ~strcmpi(surfType, a.GeometryType)
        error('Surface is not of required type: %s', surfType);
    end
    gii.DataSpace = a.DataSpace;
    if all(a.Vertices(:,3)==0), a.Vertices = a.Vertices(:,[3 1 2]); end % flat
    if     strcmp(a.AnatomicalStructurePrimary, 'CortexLeft')
        gii.Vertices{1} = a.Vertices; gii.Faces{1} = a.Faces;
    elseif strcmp(a.AnatomicalStructurePrimary, 'CortexRight')
        gii.Vertices{2} = a.Vertices; gii.Faces{2} = a.Faces;
    end
end
 
%% Return gii struct with DataSpace, Vertices and Faces.
function gii = read_gii(fname)
xml = fileread(fname); % text file basically
for i = regexp(xml, '<DataArray[\s>]') % often 2 of them
    c = regexp(xml(i:end), '.*?</DataArray>', 'match', 'once');
    [Data, i0] = gii_element(c, 'Data'); % suppose '<Data' is last element
    c = regexprep(c(1:i0-1), '<!\[CDATA\[(.*?)\]\]>', '$1'); % rm CDADA thing
    
    a = gii_attr(c, 'DataType');
    if     strcmp(a, 'NIFTI_TYPE_FLOAT32'), dType = 'single';
    elseif strcmp(a, 'NIFTI_TYPE_INT32'),   dType = 'int32';
    elseif strcmp(a, 'NIFTI_TYPE_UINT8'),   dType = 'uint8';
    else, error('Unknown GIfTI DataType: %s', a);
    end
    
    nDim = gii_attr(c, 'Dimensionality', 1);
    dim = ones(1, nDim);
    for j = 1:nDim, dim(j) = gii_attr(c, sprintf('Dim%g', j-1), 1); end
        
    Endian = gii_attr(c, 'Endian'); % LittleEndian or BigEndian
    Endian = lower(Endian(1));
 
    Encoding = gii_attr(c, 'Encoding');
    if any(strcmp(Encoding, {'Base64Binary' 'GZipBase64Binary'}))
        Data = matlab.net.base64decode(Data); % since 2016b
        if strcmp(Encoding, 'GZipBase64Binary') % HCP uses this
            Data = nii_tool('LocalFunc', 'gunzip_mem', Data);
        end
        Data = typecast(Data, dType);
        if Endian == 'b', Data = swapbytes(Data); end
    elseif strcmp(Encoding, 'ASCII') % untested
        Data = str2num(Data);
    elseif strcmp(Encoding, 'ExternalFileBinary') % untested
        nam = gii_attr(c, 'ExternalFileName');
        if isempty(fileparts(nam)), nam = fullfile(fileparts(fname), nam); end
        fid = fopen(nam, 'r', Endian);
        if fid==-1, error('ExternalFileName %s not exists'); end
        fseek(fid, gii_attr(c, 'ExternalFileOffset', 1), 'bof');
        Data = fread(fid, prod(dim), ['*' dType]);
        fclose(fid);
    else, error('Unknown Encoding: %s', Encoding);
    end
    
    if nDim>1
        if strcmp(gii_attr(c, 'ArrayIndexingOrder'), 'RowMajorOrder')
            Data = reshape(Data, dim(nDim:-1:1));
            Data = permute(Data, nDim:-1:1);
        else
            Data = reshape(Data, dim);
        end
    end
    
    Intent = gii_attr(c, 'Intent');
    if strcmp(Intent, 'NIFTI_INTENT_TRIANGLE')
        gii.Faces = Data; % 0-based
        continue;
    elseif ~strcmp(Intent, 'NIFTI_INTENT_POINTSET') % store Data only for now
        if ~isfield(gii, 'Data'), gii.Data = []; gii.Intent = []; end
        gii.Intent{end+1} = Intent;
        gii.Data{end+1} = Data;        
        continue;
    end
    
    % now only for NIFTI_INTENT_POINTSET
    meta = @(k)regexp(c, ['(?<=>' k '<.*?<Value>).*?(?=</Value>)'], 'match', 'once');
    gii.AnatomicalStructurePrimary = meta('AnatomicalStructurePrimary');
    gii.AnatomicalStructureSecondary = meta('AnatomicalStructureSecondary');
    gii.GeometricType = meta('GeometricType');
    frms = {'NIFTI_XFORM_UNKNOWN' 'NIFTI_XFORM_SCANNER_ANAT' ...
        'NIFTI_XFORM_ALIGNED_ANAT' 'NIFTI_XFORM_TALAIRACH' 'NIFTI_XFORM_MNI_152'};
    gii.DataSpace = find(strcmp(gii_element(c, 'DataSpace'), frms)) - 1;
    gii.TransformedSpace = find(strcmp(gii_element(c, 'TransformedSpace'), frms)) - 1;
    gii.MatrixData = gii_element(c, 'MatrixData', 1);
    gii.Vertices = Data;
end
 
%% Return cii/gii attribute
function val = gii_attr(ch, key, isnum)
val = regexp(ch, ['(?<=' key '=").*?(?=")'], 'match', 'once');
if nargin>2 && isnum, val = str2num(val); end %#ok<*ST2NM>
 
%% Return cii/gii element
function [val, i0] = gii_element(ch, key, isnum)
i0 = regexp(ch, ['<' key '[\s>]'], 'once');
val = regexp(ch(i0:end), ['(?<=<' key '.*?>).*?(?=</' key '>)'], 'match', 'once');
if nargin>2 && isnum, val = str2num(val); end
 
%% Open surface view or add cii to it
function cii_view(hsN)
fh = hsN.fig.UserData;
if isempty(fh) || ~ishandle(fh) % create surface figure
    fh = uifigure('HandleVisibility', 'Callback');
    cMenu = uicontextmenu(fh);
    uimenu(cMenu, 'Label', 'Reset view', 'Callback', {@cii_view_cb "reset"});
    uimenu(cMenu, 'Label', 'Zoom in' ,   'Callback', {@cii_view_cb "zoomG"});
    uimenu(cMenu, 'Label', 'Zoom out',   'Callback', {@cii_view_cb "zoomG"});
    uimenu(cMenu, 'Label', 'Change cortex color', ...
        'Callback', {@cii_view_cb "cortexColor"}, 'Separator', 'on');
    uimenu(cMenu, 'Label', 'Change surface', 'Callback', {@cii_view_cb "changeSurface"});
    saveAs = findobj(hsN.fig, 'Type', 'uimenu', 'Label', 'Save figure as');
    m = copyobj(saveAs, cMenu);
    m.Separator = 'on';
    m = m.Children; delete(m(1:3)); m = m(4:end); % pdf/eps etc too slow
    set(m, 'Callback', {@nii_viewer_cb "save" fh});
    if ispc || ismac
        uimenu(cMenu, 'Label', 'Copy figure', 'Callback', {@nii_viewer_cb 'copy' fh});
    end
    
    r = 0.96; % width of two columns, remaining for colorbar
    pos = [0 1 r 1; 0 0 r 1; r 1 r 1; r 0 r 1] / 2;
    gii = cii2nii(); % get buffered Anatomical gii
    hs.frame = uipanel(fh, 'Units', 'normalized', 'Position', [0 0 1 1], ...
        'BackgroundColor', 'k', 'UIContextMenu', cMenu);
    for ig = 1:2
        v = gii.Vertices{ig};
        lim = [min(v)' max(v)'];
        im = ones(size(v,1), 3, 'single') * 0.667;
        for i = ig*2+[-1 0]
            hs.ax(i) = axes(hs.frame, 'Position', pos(i,:), 'CameraViewAngle', 6.8);
            axis vis3d; axis equal; axis off;
            set(hs.ax(i), 'XLim', lim(1,:), 'YLim', lim(2,:), 'ZLim', lim(3,:));
            hs.ax(i).Toolbar.Visible = 'off';
            hs.patch(i) = patch(hs.ax(i), 'EdgeColor', 'none', ...
                'Faces', gii.Faces{ig}+1, 'Vertices', v, ...
                'FaceVertexCData', im, 'FaceColor', 'interp', 'FaceLighting', 'gouraud');
            hs.light(i) = camlight('infinite'); material dull;
        end
    end
    set(hs.patch, 'ButtonDownFcn', {@cii_view_cb "buttonDownPatch"}, 'UIContextMenu', cMenu);
    
    hs.ax(5) = axes(hs.frame, 'Position', [r 0.1 1-r 0.8], 'Visible', 'off');
    hs.colorbar = colorbar(hs.ax(5), 'PickableParts', 'none', ...
        'Location', 'East', 'Visible', hsN.colorbar.Visible);
            
    fh.Position(3:4) = [1/r diff(lim(3,:))/diff(lim(2,:))] * 600 + 4;
    srn = get(0, 'MonitorPositions');
    dz = srn(1,4)- 60 - sum(fh.Position([2 4]));
    if dz<0, fh.Position(2) = fh.Position(2) + dz; end
 
    fh.WindowButtonUpFcn   = {@cii_view_cb "buttonUp"};
    fh.WindowButtonDownFcn = {@cii_view_cb "buttonDown"};
    fh.WindowButtonMotionFcn={@cii_view_cb "buttonMotion"};
    fh.WindowKeyPressFcn   = {@cii_view_cb "keyFcn"};
    fh.UserData = struct('xy', [], 'xyz', [], 'hemi', [], 'deg', [-90 0], ...
                         'color', [1 1 1]*0.667);
    hsN.fig.UserData = fh;
    hs.gii = gii.Vertices;
    hs.fig = fh;
    hs.hsN = hsN;
    guidata(fh, hs);
    set_cii_view(hs, [-90 0]);
    set(hs.patch, 'VertexNormalsMode', 'auto');
    drawnow; set(hs.patch, 'VertexNormalsMode', 'manual'); % faster update
end
 
set_colorbar(hsN);
cii_view_cb(fh, [], "volume");
cii_view_cb(fh, [], "background");
 
%% cii_view callbacks 
function cii_view_cb(h, ev, cmd)
if isempty(h) || ~ishandle(h), return; end
hs = guidata(h);
if cmd == "buttonMotion" % Rotate gii and set light direction
    if isempty(hs.fig.UserData.xy), return; end % button not down
    d = hs.fig.CurrentPoint - hs.fig.UserData.xy; % location change
    d = hs.fig.UserData.deg - d; % new azimuth and elevation
    set_cii_view(hs, d);
    hs.fig.UserData.xyz = []; % so not set cursor in nii_viewer
elseif cmd == "buttonDownPatch" % get button xyz
    if ~strcmpi(hs.fig.SelectionType, 'normal') % button 1 only
        hs.fig.UserData.xyz = []; return;
    end
    hs.fig.UserData.hemi = 1 + any(h.Parent == hs.ax(3:4));
    hs.fig.UserData.xyz = ev.IntersectionPoint;
elseif cmd == "buttonDown" % figure button down: prepare for rotation
    if ~strcmpi(hs.fig.SelectionType, 'normal'), return; end % not button 1
    hs.fig.UserData.xy = hs.fig.CurrentPoint;
    hs.fig.UserData.deg = hs.ax(1).View;
elseif cmd == "buttonUp" % set cursor in nii_viewer
    hs.fig.UserData.xy = []; % disable buttonMotion
    xyz = hs.fig.UserData.xyz;
    if isempty(xyz), return; end
    ig = hs.fig.UserData.hemi;
    v = hs.patch(ig*2).Vertices;
    if ~isequal(v, hs.gii{ig})
        [~, iv] = min(sum((v - xyz) .^ 2, 2));
        xyz = hs.gii{ig}(iv,:);
    end
    c = round(hs.hsN.bg.Ri * [double(xyz(:)); 1]) + 1;
    for i = 1:3, hs.hsN.ijk(i).Value = c(i); end
    nii_viewer_cb([], [], 'ijk', hs.hsN.fig);
elseif cmd == "reset"
    set_cii_view(hs, [-90 0]);
    set(hs.ax(1:4), 'CameraViewAngle', 6.8);
elseif cmd == "cortexColor"
    c = uisetcolor(hs.fig.UserData.color, 'Set cortical surface color');
    if numel(c) ~= 3, return; end
    hs.fig.UserData.color = c;
    set_cdata_cii(hs);
elseif cmd == "changeSurface"
    gii = get_surfaces(size(hs.gii{1}, 1));
    if isempty(gii), return; end
    if size(gii.Vertices{1},1) ~= size(hs.gii{1},1)
        errordlg('GIfTI has different number of vertices from CIfTI');
        return;
    end
    ind = [isequal(gii.Vertices{1}, hs.patch(1).Vertices) ...
           isequal(gii.Vertices{2}, hs.patch(3).Vertices)];
    if all(ind), return; end % no change
    
    for i = find(~ind)
        ip = i*2+(-1:0);
        set(hs.patch(ip), 'Vertices', gii.Vertices{i}, 'Faces', gii.Faces{i}+1);
        lim = [min(gii.Vertices{i})' max(gii.Vertices{i})'];
        set(hs.ax(ip), 'ZLim', lim(3,:), 'YLim', lim(2,:));
        try set(hs.ax(ip), 'XLim', lim(1,:)); end % avoid error for flat
    end
    
    set_cdata_cii(hs);
    try
        set(hs.patch, 'VertexNormalsMode', 'auto');
        drawnow; set(hs.patch, 'VertexNormalsMode', 'manual');
    end
elseif cmd == "zoomG"
    va = hs.ax(1).CameraViewAngle;
    if strcmp(h.Label, 'Zoom in'), va = va * 0.9;
    else, va = va * 1.1;
    end
    set(hs.ax(1:4), 'CameraViewAngle', va);
elseif cmd == "keyFcn"
    if any(strcmp(ev.Key, ev.Modifier)), return; end % only modifier
    figure(hs.fig); % avoid focus to Command Window
    if ~isempty(intersect({'control' 'command'}, ev.Modifier))
        if any(strcmp(ev.Key, {'add' 'equal'}))
            h = findobj(hs.fig, 'Type', 'uimenu', 'Label', 'Zoom in');
            cii_view_cb(h, [], 'zoomG');
        elseif any(strcmp(ev.Key, {'subtract' 'hyphen'}))
            h = findobj(hs.fig, 'Type', 'uimenu', 'Label', 'Zoom out');
            cii_view_cb(h, [], 'zoomG');
        elseif any(strcmp(ev.Key, {'a' 'r'}))
            figure(hs.hsN.fig);
            key = java.awt.event.KeyEvent.(['VK_' upper(ev.Key)]);
            java.awt.Robot().keyPress(key);
            java.awt.Robot().keyRelease(key);
        end
    elseif any(strcmp(ev.Key, {'space' 'comma' 'period'}))
        KeyPressFcn(hs.hsN.fig, ev);
    end
    
% Rest not called by surface callback, but from nii_viewer_cb
elseif any(cmd == ["lb" "ub" "lut" "toggle" "alpha" "volume" "stack" "close" "closeAll"])
    set_cdata_cii(hs);
    if cmd == "volume", cii_view_cb(h, [], "file"); end
elseif cmd == "file" % update MapName if available
    p = get_para(hs.hsN);
    try mName = p.nii.NamedMap{p.volume}.MapName; catch, mName = ''; end
    frm = formcode2str(p.nii.hdr.sform_code);
    hs.fig.Name = ['cii_view - ' mName ' (' frm ')'];
elseif cmd == "background"
    clr = hs.hsN.frame.BackgroundColor;
    hs.frame.BackgroundColor = clr;
    hs.colorbar.EdgeColor = 1-clr;
elseif cmd == "colorbar" % colorbar on/off
    hs.colorbar.Visible = hs.hsN.colorbar.Visible;
end
 
%% set surface FaceVertexCData
function set_cdata_cii(hs)
imP{1} = ones(size(hs.gii{1},1), 1, 'single') * hs.fig.UserData.color;
imP{2} = ones(size(hs.gii{2},1), 1, 'single') * hs.fig.UserData.color;
for j = size(hs.hsN.files.Data,1):-1:1
    p = get_para(hs.hsN, j);
    if ~p.show || ~isfield(p.nii, 'cii'), continue; end
    for i = 1:2 % hemispheres
        [im, alfa] = lut2img(p.nii.cii{i}(:, p.volume), p);
        alfa = p.alpha * single(alfa>0);
        im = permute(im, [1 3 2]); % nVertices x 3
        imP{i} = imP{i} .* (1-alfa) + im .* alfa;
    end
end
set(hs.patch(1:2), 'FaceVertexCData', imP{1});
set(hs.patch(3:4), 'FaceVertexCData', imP{2});
% tic; drawnow; toc
 
%% set surface view
function set_cii_view(hs, ae)
ae = [ae; ae(1)+180 -ae(2); -ae(1) ae(2); -ae(1)-180 -ae(2)];
for i = 1:4
    hs.ax(i).View  = ae(i,:);
    camlight(hs.light(i), 'headlight');
end
 
%% normalize columns
function v = normc(M)
v = M ./ sqrt(sum(M .* M));
 
%% reorient nii to diagnal major
function [nii, perm, flp] = nii_reorient(nii, leftHand, ask_code)
if nargin<3, ask_code = []; end
[R, frm] = nii_xform_mat(nii.hdr, ask_code);
dim = nii.hdr.dim(2:4);
pixdim = nii.hdr.pixdim(2:4);
[R, perm, flp] = reorient(R, dim, leftHand);
fps = bitand(nii.hdr.dim_info, [3 12 48]) ./ [1 4 16];
if ~isequal(perm, 1:3)
    nii.hdr.dim(2:4) = dim(perm);
    nii.hdr.pixdim(2:4) = pixdim(perm);
    nii.hdr.dim_info = [1 4 16] * fps(perm)' + bitand(nii.hdr.dim_info, 192);
    nii.img = permute(nii.img, [perm 4:8]);
end
sc = nii.hdr.slice_code;
if sc>0 && flp(fps==3)
    nii.hdr.slice_code = sc+mod(sc,2)*2-1; % 1<->2, 3<->4, 5<->6
end
if isequal(perm, 1:3) && ~any(flp), return; end
if frm(1) == nii.hdr.sform_code % only update matching form
    nii.hdr.sform_mat = R(1:3,:);
end
if frm(1) == nii.hdr.qform_code
    nii.hdr.qform_mat = R(1:3,:);
    nii.hdr.qoffset_xyz = R(1:3,4);
    R0 = normc(R(1:3, 1:3));
    dcm2quat = dicm2nii('', 'dcm2quat', 'func_handle');
    [q, nii.hdr.pixdim(1)] = dcm2quat(R0);
    nii.hdr.quatern_bcd = q(2:4);
end
for i = find(flp), nii.img = flip(nii.img, i); end
 
%% Change stack order for uifigure
function im_stack(ax, ind)
typ = arrayfun(@(a)a.Type, ax(1).Children, 'UniformOutput', false);
iImg = find(typ=="image" | typ=="quiver");
for j = 1:3
    ch = ax(j).Children;
    ax(j).Children(iImg) = ch(iImg(ind));
end
 
%% Robot mouse to remove focus from a component
function focus_off(fh)
if ~isequal(ancestor(gcbo, 'figure'), fh), return; end % call from cii_view
xy = get(0, 'PointerLocation'); % for later restore
pos = getpixelposition(fh);
tr = pos(1:2) + pos(3:4) - 4; % close to top-right corner
try
    rob = java.awt.Robot();
    if ismac
        res = get(0, 'ScreenSize');
        rob.mouseMove(tr(1), res(4)-tr(2));
        rob.mousePress(16); rob.mouseRelease(16); % 16=BUTTON1
        rob.mouseMove(xy(1), res(4)-xy(2));
    else
        set(0, 'PointerLocation', tr);
        rob.mousePress(16); rob.mouseRelease(16);
        set(0, 'PointerLocation', xy);
    end
end
 
%% update uitable width
function table_width(h)
x = max(cellfun(@numel, h.Data(:,2))) * 6; % pixels per character
if size(h.Data,1)>3, x = x + 18; end % vertical scroll bar
ph = ancestor(h, 'uipanel');
x = max(60, min(x+48, ph.Position(3)-474)); % 474 = min param panel
h.Parent.ColumnWidth{1} = x; % uigridlayout {x '1x'}
x2 = x - 24; if size(h.Data,1)>3, x2 = x2 - 18; end
h.ColumnWidth{2} = x2;

%% Return true if input is char or single string (R2016b+)
function tf = ischar(A)
tf = builtin('ischar', A) || isStringScalar(A);

%% manual alignment GUI
function d = manual_align(h, ~)
fh = ancestor(h, 'figure');
pos = getpixelposition(fh);
d = uifigure('Name', 'Manual alignment', 'HandleVisibility', 'callback', ...
    'Position', [pos(1)+pos(3) pos(2) 334 210], 'Resize', 'off', 'UserData', fh);
uitxt = @(p3,str)uilabel(d, 'Text', str, 'Position', [p3-[0 4 0] 22], ...
    'Tooltip', str, 'HorizontalAlignment', 'right');
spinner = @(x,y,v,lim,step)uispinner(d, 'Position', [x y 60 24], 'Value', v, ...
    'ValueDisplayFormat', '%.3g', 'Limits', lim, 'Step', step, 'ValueChangedFcn', @align_cb);

x = [0 72 144] + 120;
y = (0:30:120) + 60;
uitxt([x(1) y(5) 28], 'X'); uitxt([x(2) y(5) 28], 'Y'); uitxt([x(3) y(5) 28], 'Z');
uitxt([12 y(4) 98], 'Translation (mm)');
for i = 1:3, hs.p(i) = spinner(x(i), y(4), 0, [-99 99], 1); end
uitxt([12 y(3) 98], 'Rotation (deg)');
for i = 4:6, hs.p(i) = spinner(x(i-3), y(3), 0, [-45 45], 0.5); end
uitxt([12 y(2) 98], 'Scaling');
for i = 7:9, hs.p(i) = spinner(x(i-6), y(2), 1, [0.5 2], 0.01); end
uitxt([12 y(1) 98], 'Shearing');
for i = 10:12, hs.p(i) = spinner(x(i-9), y(1), 0, [-1 1], 0.01); end

uibtn = @(x,w,btn)uibutton(d, 'Text', btn, 'Position', [x 18 w 24], 'ButtonPushedFcn', @align_cb);
uibtn(18, 52, 'Reset'); uibtn(102, 90, 'vector 1x12'); uibtn(224, 90, 'matrix 4x4');
guidata(d, hs.p);

%% manual alignment callback
function align_cb(h, ~)
p = guidata(h);
v12 = [p.Value];
if h.Type == "uibutton"
    if h.Text == "Reset"
        v12 = [0 0 0  0 0 0  1 1 1  0 0 0];
        for i = 1:12, p(i).Value= v12(i); end
    elseif h.Text == "vector 1x12"
        fprintf('v12 = ['); fprintf('%.3g ', v12); fprintf('\b];\n');
        return;
    elseif h.Text == "matrix 4x4"
        M = affine_mat(v12);
        global MANUAL_ALIGN; MANUAL_ALIGN = M; %#ok
        fprintf('m44 = [\n'); fprintf('%9.3g%9.3g%9.3g%9.3g\n', M'); fprintf('];\n');
        return;
    end
end
hs = guidata(get(ancestor(h, 'figure'), 'UserData'));
p = get_para(hs);
p.hsI(1).UserData.Ri = inv(affine_mat(v12) * p.hsI(1).UserData.R);
set_cdata(hs);
 
%% Return transformation matrix from a vector of 3, 6, 9 or 12
function A = affine_mat(P)
% Return 4x4 affine transformation matrix
% P(1:3)   - xyz translation
% P(4:6)   - xyz rotation about pitch, roll, yaw (degrees)
% P(7)     - single scaling for xyz if numel(P)==7
% P(7:9)   - xyz scaling if numel(P)>=9
% P(10:12) - xyz shearing

A = eye(4); A(1:3, 4) = P(1:3);
if numel(P)==3, return; end
ca = cosd(P(4:6)); sa = sind(P(4:6));
A(1:3, 1:3) = [1 0 0; 0 ca(1) sa(1); 0 -sa(1) ca(1)] * ...
              [ca(2) 0 sa(2); 0 1 0; -sa(2) 0 ca(2)] * ...
              [ca(3) sa(3) 0; -sa(3) ca(3) 0; 0 0 1];
if numel(P)==6, return; end
if numel(P)==7, A(1:3, 1:3) = A(1:3, 1:3) * P(7); return; end
A(1:3, 1:3) = A(1:3, 1:3) * diag(P(7:9));
if numel(P)==9, return; end
S = eye(4); S([5 9 10]) = P(10:12);
A = A * S;

%%
