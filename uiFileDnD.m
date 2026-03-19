function uiFileDnD(obj, dropFcn)
% Set up a callback when file/folder is dropped onto a (ui)figure.
% 
% The obj can be figure/uifigure or its component.
% 
% dropFcn is the callback function when a file is dropped. Its syntax is the
% same as general Matlab callback, like @myFunc or {@myFunc myOtherInput}.
% In the callback, the first argument is the obj, and 2nd a struct containing
%      dropXY: [x y] % drop location in pixels relative to figure
%     ctrlKey: 0 % true if Ctrl key is down while dropping
%    shiftKey: 0 % true if Shift key is down while dropping
%       names: {'/myPath/myFile'} % cellstr for full file/folder names
% 
% Example to show dropped file/folder onto uilistbox of uifigure:
%  obj = uilistbox(uifigure, 'Position', [80 100 400 100]);
%  uiFileDnD(obj, @(o,dat)set(o,'Items',dat.names));
%
% Example to show dropped file/folder onto listbox of figure:
%  obj = uicontrol(figure, 'Style', 'listbox', 'Position', [80 100 400 100]);
%  uiFileDnD(obj, @(o,dat)set(o,'String',dat.names));

% 201001 Wrote it, by Xiangrui.Li at gmail.com 
% 201023 Remove uihtml by using ww.executeJS
% 260125 need ForceIndependentlyHostedFigures for R2025+. Thx EricMagalhaesDelgado 
% 260129 Rename from DnD_uifigure since works for figure too
% 260305 Add 2nd hidden button to update rects, so work for text drop too.

narginchk(2, 2);
if isempty(obj), obj = uifigure; end
if numel(obj)>1 || ~ishandle(obj)
    error('uiFileDnD:badInput', 'obj must be a single (ui)figure component');
end

fh = ancestor(obj, 'figure');
old = warning('off'); resetWarn = onCleanup(@()warning(old)); % MATLAB:structOnObject
drawnow; fhS = struct(fh);

% This if-end block is for figure() before R2025a and can be removed in the future,
% together with java_dnd.m & MLDropTarget.class
if ~isfield(fhS, 'Controller') || isempty(fhS.Controller)
    if ~exist("java_dnd.m", "file")
        error("You are using figure() under R2024b or earlier. Get the package at:" + ...
            newline + "https://github.com/xiangruili/uiFileDnD");
    end
    java_dnd(obj, dropFcn);
    return
end

hBtn = findall(fh, 'Type', 'uibutton', 'Tag', 'uiFileDropBtn');
if ~isempty(hBtn)
    ind = find([hBtn.UserData{:,2}]==obj, 1, 'last');
    if isempty(ind), hBtn.UserData(end+1,:) = {dropFcn obj};
    else, hBtn.UserData{ind,1} = dropFcn;
    end
    return;
end

try
    ww = struct(struct(fhS.Controller).PlatformHost).CEF;
    ww.enableDragAndDropAll; % DnD to whole figure: no-op for Linux till 2024b
catch me
    if regexp(me.message, 'Unrecognized.+CEF') % >=R2025a
        error("For Matlab R2025a or later, add the following line into " + ...
            "your startup.m file (create it if not exists):" + newline +...
            "try addprop(groot, 'ForceIndependentlyHostedFigures'); catch, end");
    elseif regexp(me.message, 'Unrecognized.+enableDragAndDropAll') % <R2020b
        error('Matlab R2020b or later needed for file drag and drop onto uifigure');
    else % other error
        rethrow(me);
    end
end
hBtn = uibutton(fh, 'Position', [1 1 0 0], 'Text', 'uiFileDropBtn', ...
    'ButtonPushedFcn', {@drop ww}, 'Visible', 'off', 'HandleVisibility', 'off');
h = copyobj(hBtn, fh); 
set(h, 'Text', 'uiFileDragBtn', 'ButtonPushedFcn', {@dragEnter ww hBtn});
set(hBtn, 'Tag', 'uiFileDropBtn', 'UserData', {dropFcn obj});
jsStr = fileread(mfilename("fullpath")+".m");
jsStr = regexp(jsStr, 'javascript\s+\%\{\s+(.*?)\%\}', 'tokens', 'once');
drawnow; ww.executeJS(jsStr{1});
ww.FileDragDropCallback = @(o,names)set(hBtn,'Text',cellstr(names));

%% fired by javascript fake button press
function dragEnter(~, ~, ww, h)
h.UserData = h.UserData(isvalid([h.UserData{:,2}]), :); % remove invalid
for i = size(h.UserData,1):-1:1 % redo in case pos changed or resized
    obj = h.UserData{i,2};
    p{i} = round(getpixelposition(obj, true));
    if p{i}(4)==0, p{i} = getpixelposition(obj.Parent, true); % axes' child
    elseif obj.Type == "figure", p{i}(1:2) = 1;
    end
end
ww.executeJS(['rects=' jsonencode(p)]);

function drop(h, ~, ww)
dat = jsondecode(ww.executeJS('Data'));
dat.dropXY = dat.dropXY';
if isempty(dat.names), dat.names = h.Text; end
args = [h.UserData(dat.index+1,:) rmfield(dat, 'index')];
if iscell(args{1}), args = [args{1}(1) args(2:3) args{1}(2:end)]; end
feval(args{:});

%% following script must start with % javascript, will run once for 1st obj
% javascript
%{
let rects = [];
let Data = {dropXY: [0,0], ctrlKey: false, shiftKey: false, names: [], index: 0};
const bDrag = [...document.querySelectorAll('.mwPushButton')].find(btn => btn.textContent.trim() === 'uiFileDragBtn');
const bDrop = [...document.querySelectorAll('.mwPushButton')].find(btn => btn.textContent.trim() === 'uiFileDropBtn');

document.ondragenter = (e) => {
    bDrag.click(); // ask to update rects
    if (e.dataTransfer) e.dataTransfer.dropEffect = "none";
    return false;
};
document.ondragover = (e) => {
    e.returnValue = false;
    if (!e.dataTransfer) return;
    const x = e.clientX+1, y = window.innerHeight-e.clientY;
    for (let i = rects.length-1; i >= 0; i--) {
        const p = rects[i]; // [left bottom width height]
        if (x>p[0] && y>p[1] && x<p[0]+p[2] && y<p[1]+p[3]) {
            Data.index = i;
            Data.dropXY = [x,y];
            return;
        };
    };
    e.dataTransfer.dropEffect = 'none'; // disable drop
};
document.ondrop = (e) => {
    e.returnValue = false;
    Data.ctrlKey = e.ctrlKey;
    Data.shiftKey = e.shiftKey;
    Data.names = [];
    if (e.dataTransfer.types.includes("text/plain")) {
       Data.names = [e.dataTransfer.getData("text/plain")];
    }
    bDrop.click();
};
%}
