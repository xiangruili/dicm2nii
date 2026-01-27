function uiFileDnD(target, dropFcn)
% Set up a callback when file/folder is dropped onto a figure/uifigure component.
% 
% The target can be figure/uifigure or its component.
% 
% dropFcn is the callback function when a file is dropped. Its syntax is the
% same as general Matlab callback, like @myFunc or {@myFunc myOtherInput}.
% In the callback, the first argument is the target, and 2nd the data containing
%     ctrlKey: 0 % true if Ctrl key is down while dropping
%    shiftKey: 0 % true if Shift key is down while dropping
%       names: {'/myPath/myFile'} % cellstr for full file/folder names
% 
% Example to show dropped file/folder onto uilistbox:
%  target = uilistbox(uifigure, 'Position', [80 100 400 100]);
%  uiFileDnD(target, @(o,dat)set(o,'Items',dat.names));

% 201001 Wrote it, by Xiangrui.Li at gmail.com 
% 201023 Remove uihtml by using ww.executeJS
% 260125 use ForceIndependentlyHostedFigures for R2025+. Thx EricMagalhaesDelgado 
% 260126 Rename from DnD_uifigure since works for figure too

narginchk(2, 2);
if isempty(target), target = uifigure; end
if numel(target)>1 || ~ishandle(target)
    error('uiFileDnD:badInput', 'target must be a single (ui)figure component');
end

fh = ancestor(target, 'figure');

% this block and java_dnd.m&MLDropTarget.class can be removed in future
try isUIFigure = fh.isUIFigure; catch, isUIFigure = false; end
if ~isUIFigure && verLessThan('matlab', '25.1') %#ok
    java_dnd(target, dropFcn);
    return
end

hBtn = findall(fh, 'Type', 'uibutton', 'Tag', 'uiFileDnD');
if ~isempty(hBtn)
    hBtn.UserData(end+1,:) = {dropFcn target};
    return;
end

drawnow;
old = warning('off'); % MATLAB:structOnObject
cln = onCleanup(@()warning(old));
try
    ww = struct(struct(struct(fh).Controller).PlatformHost).CEF;
    ww.enableDragAndDropAll; % DnD to whole uifigure: no-op for Linux till 2021a
catch me
    if verLessThan('matlab', '9.9') %#ok < R2020b
        error('Matlab R2020b or later needed for file drag and drop');
    elseif verLessThan('matlab', '25.1') %#ok < R2025a 
        rethrow(me);
    else
        error("For Matlab R2025a or later, add the following line into " + ...
            "your startup.m file (create it if not exists):" + newline +...
            "try addprop(groot, 'ForceIndependentlyHostedFigures'); catch, end");
    end
end
hBtn = uibutton(fh, 'Position', [1 1 0 0], 'Text', '4JS2identify_me', ...
    'ButtonPushedFcn', {@drop ww}, 'UserData', {dropFcn target}, ...
    'Tag', 'uiFileDnD', 'Visible', 'off', 'HandleVisibility', 'off');

jsStr = char(strjoin([ ... % webwindow accepts only char at least for R2020b
    % """use strict"";"
    "let uiFileDnD = {rects: [], lastOver: 0,"
    "   data: {ctrlKey: false, shiftKey: false, index: 0},"
    "   button: [...document.querySelectorAll('.mwPushButton')].find("
    "      btn => btn.textContent.trim() === '"+hBtn.Text+"')};"
    "document.ondragenter = (e) => { // prevent default before firing ondragover"
    "  e.dataTransfer.dropEffect = 'none';"
    "  return false;"
    "};"
    "document.ondragover = (e) => {"
    "  e.returnValue = false; // preventDefault & stopPropagation"
    "  let now = new Date().getTime();"
    "  if (now < uiFileDnD.lastOver+16) { return; }"
    "  uiFileDnD.lastOver = now;"
    "  let x = e.clientX+1, y = document.body.clientHeight-e.clientY;"
    "  for (let i = uiFileDnD.rects.length-1; i >= 0; i--) {"
    "    let p = uiFileDnD.rects[i]; // [left bottom width height]"
    "    if (x>=p[0] && y>=p[1] && x<p[0]+p[2] && y<p[1]+p[3]) {"
    "      uiFileDnD.data.index = i; // target index in rects"
    "      return; // keep OS default dropEffect"
    "    };"
    "  };"
    "  e.dataTransfer.dropEffect = 'none'; // disable drop"
    "};"
    "document.ondrop = (e) => {"
    "  e.returnValue = false;"
    "  uiFileDnD.data.ctrlKey = e.ctrlKey;"
    "  uiFileDnD.data.shiftKey = e.shiftKey;"
    "  uiFileDnD.button.click(); // fire Matlab callback"
    "};" ], newline));
drawnow; ww.executeJS(jsStr);
ww.FileDragDropCallback = {@dragEnter hBtn};

%% fired when drag enters uifigure
function dragEnter(ww, names, h)
for i = size(h.UserData,1):-1:1 % redo in case pos changed or resized
    p{i} = round(getpixelposition(h.UserData{i,2}, 1));
    if h.UserData{i,2}.Type == "figure", p{i}(1:2) = 1; end
end
ww.executeJS(['uiFileDnD.rects=' jsonencode(p)]);
h.Text = cellstr(names); % store file names

%% fired by javascript fake button press in ondrop
function drop(h, ~, ww)
dat = jsondecode(ww.executeJS('uiFileDnD.data'));
dat.names = h.Text;
args = [h.UserData(dat.index+1,:) rmfield(dat, 'index')];
if iscell(args{1}), args = [args{1}(1) args(2:3) args{1}(2:end)]; end
feval(args{:});
%%

%%
