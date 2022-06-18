function varargout = esVidTrack(varargin)
% ESVIDTRACK MATLAB code for esVidTrack.fig
%      ESVIDTRACK, by itself, creates a new ESVIDTRACK or raises the existing
%      singleton*.
%
%      H = ESVIDTRACK returns the handle to a new ESVIDTRACK or the handle to
%      the existing singleton*.
%
%      ESVIDTRACK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ESVIDTRACK.M with the given input arguments.
%
%      ESVIDTRACK('Property','Value',...) creates a new ESVIDTRACK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before esVidTrack_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to esVidTrack_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help esVidTrack

% Last Modified by GUIDE v2.5 26-Mar-2014 13:14:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @esVidTrack_OpeningFcn, ...
    'gui_OutputFcn',  @esVidTrack_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before esVidTrack is made visible.
function esVidTrack_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to esVidTrack (see VARARGIN)

% Choose default command line output for esVidTrack
handles.output = hObject;

handles.fn_path = mfilename('fullpath');
if exist (fullfile(fileparts(handles.fn_path),'esVidTrack_data.mat'),'file')
    load(fullfile(fileparts(handles.fn_path),'esVidTrack_data.mat'));
    info = get(mds,'UserData');
    mask_arena = poly2mask(info.position.arena(:,1), info.position.arena(:,2), info.Height, info.Width);
    setappdata(handles.figure1,'mask_arena',mask_arena);
    handles.position = info.position;
else
    mds = dataset;
end

setappdata(handles.figure1, 'mds', mds); % master dataset

handles.skip_frame = str2num(get(handles.edit_skip_frame,'String'));
handles.level = str2num(get(handles.edit_level,'String'));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes esVidTrack wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = esVidTrack_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in push_load.
function push_load_Callback(hObject, eventdata, handles)
% hObject    handle to push_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile({'*.*;*','All files'});

vobj = VideoReader(fullfile(pathname, filename));
handles.info = parse_filename(filename); % initialize info
handles.info.pathname = pathname;

mds = getappdata(handles.figure1, 'mds');

if ~isempty(mds)
    filenames = unique(mds.filename);
    if any(strcmp(filenames,handles.info.filename ))
        uiwait(msgbox('Warning. File already tracked'))
    end
end

handles.info.numFrames = get(vobj,'NumberOfFrames');
handles.info.frameRate = vobj.FrameRate;
handles.info.Height = vobj.Height;
handles.info.Width = vobj.Width;
handles.vindex = 1;

img1 = read(vobj,handles.vindex);

if isappdata(handles.figure1, 'mask_arena')
    mask_arena = getappdata(handles.figure1, 'mask_arena');
    imshowpair(img1,~mask_arena)
else
    imshow(img1);
end

tdata = [fieldnames(handles.info) struct2cell(handles.info)];
set(handles.uitable1,'Data',tdata)

handles.vobj = vobj;
guidata(hObject, handles);

function info = parse_filename(filename)
% important: Filename must look like this: 031014 WT1 W-W T0.mov
% only length of W-W can vary

[~,fname,ext] = fileparts(filename);
spaces = regexp(filename,'\s');
dash = regexp(filename,'-');

info.filename = fname;
info.ext = ext;
info.date = fname(1:spaces(1)-1);
info.animal = fname(spaces(1)+1:spaces(2)-1);
info.exp_left = fname(spaces(2)+1:dash-1);
info.exp_right = fname(dash+1:spaces(3)-1);
info.time = fname(spaces(3)+1:end);


% --- Executes on button press in push_set_arena.
function push_set_arena_Callback(hObject, eventdata, handles)
% hObject    handle to push_set_arena (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);
himg = imshow(read(handles.vobj,1));
% uiwait(msgbox('Outline Arena'))

display('Outline Arena')
he_arena = impoly(gca); % create a draggable, resizable polygon
handles.position.arena = wait(he_arena); % wait for a double click on polygon
mask_arena = createMask(he_arena,himg); % create mask from polygon
setappdata(handles.figure1, 'mask_arena', mask_arena)

delete(he_arena)
imshowpair(read(handles.vobj,handles.vindex),...
    ~mask_arena);


% uiwait(msgbox('Draw line to Divide Arena'))
display('Draw line to Divide Arena')
he_line = imline(gca); % create a draggable, resizable polygon
handles.position.line = wait(he_line); % wait for a double click on polygon
delete(he_line);


hold on
plot(handles.position.line(:,1), handles.position.line(:,2),'w');

guidata(hObject, handles);

% --- Executes on button press in push_track.
function push_track_Callback(hObject, eventdata, handles)
% hObject    handle to push_track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mds = getappdata(handles.figure1, 'mds');

if ~isempty(mds) && any(strcmp(mds.filename, handles.info.filename ))
    button = questdlg('Same filename. Overwrite previously captured data?','Warning','yes','no','no');
    switch button
        case 'yes'
            mds(strcmp(mds.filename,handles.info.filename ),:) = [];
        case 'no'
            return
    end
end

count = handles.info.numFrames;
skip_frame = handles.skip_frame;
comp_frame = floor(skip_frame/2);

frame_idx = 1:skip_frame:count-skip_frame-1;
if isappdata(handles.figure1,'mask_arena')
    mask_arena = getappdata(handles.figure1, 'mask_arena');
else
    beep
    display ('set mask arena first')
    return
end

hw = waitbar(0,'Reading frames');
mouse_centroid = zeros(numel(frame_idx), 2);
axes(handles.axes1);
level = handles.level; 
for n = 1:numel(frame_idx)
    
    % read images
    img1 = read(handles.vobj,frame_idx(n));
    img2 = read(handles.vobj,frame_idx(n)+comp_frame);
    
    % set image 1 to grayscale and complement   
    imgg1 = single(mat2gray(rgb2gray(img1)));
    imgc1 = imcomplement(imgg1);
    imgc1(~mask_arena) = 0;
    
    % imgt1 = imtophat(imgc1, strel('disk',75));
    %     level = graythresh(imgc1);
    
    % threshold
    bw1 = im2bw(imgc1, level);
    bw1 = bwareaopen(bw1, 200);
    bw1 = bw1 & mask_arena;
    
    % repeat
    imgg2 = single(mat2gray(rgb2gray(img2)));
    imgc2 = imcomplement(imgg2);
    imgc2(~mask_arena) = 0;
    
    % imgt2 = imtophat(imgc2, strel('disk',75));
    %     level = graythresh(imgc2);
    bw2 = im2bw(imgc2, level);
    bw2 = bwareaopen(bw2, 200);
    bw2 = bw2 & mask_arena;
    
    % compare the complements
    imgc1(~bw1) = 0;
    imgc2(~bw2) = 0;
    imgd = imsubtract(imgc1, imgc2);
    
    % reconstruct based on differences
    bwrc = imreconstruct(logical(imgd), bw2);
%     L = bwlabel(bwrc);
    stats = regionprops(logical(bwrc),'basic'); % get centroid of all white pixels           
    [~,largest_idx] = max([stats.Area]);
    cla
    imshowpair(img2, bwrc)

    mouse_centroid(n,:) = stats(largest_idx).Centroid;
    hold on
    plot(mouse_centroid(1:n,1), mouse_centroid(1:n,2),'w')
    try
        waitbar(n/numel(frame_idx))
    catch
        display('tracking canceled')
        return
    end
end

delete(hw);
display('tracking complete')

% set right / left data points
line_position = handles.position.line;
vq = interp1(line_position(:,2),line_position(:,1), mouse_centroid(:,2)); % interpolate Y positions into an X position on the line
in_right = mouse_centroid(:,1) > vq; % logical array of whether position is in right arena

% create current dataset and save
cds = dataset({frame_idx', 'frame'}, {mouse_centroid, 'XY'}, in_right);
Info = handles.info;
Info.position = handles.position;
Info.skip_frame = skip_frame;
cds = set(cds, 'UserData',Info);
save(fullfile(handles.info.pathname , handles.info.filename ), 'cds');

total_time = length(cds) * (skip_frame / handles.info.frameRate);
time_right = sum(cds.in_right) * (skip_frame / handles.info.frameRate);
time_left = total_time - time_right;
assignin('base','cds',cds);

% update master dataset, load into workspace, and save
mds = vertcat(mds, horzcat(struct2cell2dataset(handles.info), dataset(skip_frame, time_left, time_right, total_time)));
mds = set(mds,'UserData',Info);
setappdata(handles.figure1, 'mds', mds);
save(fullfile(fileparts(handles.fn_path),'esVidTrack_data.mat'), 'mds');
assignin('base','mds',mds);

% show scatter plot
imshow(img1)
hold on
plot(mouse_centroid(:,1), mouse_centroid(:,2),'w')
gscatter(mouse_centroid(:,1), mouse_centroid(:,2), in_right)

function centroid = compare_engine(img1, img2,level)





function centroid = compare_engine_old(img1, img2, mask_arena)
% too susceptible to image flickering
img_diff = mat2gray(rgb2gray(imsubtract(img1, img2))); % subtract each channel first before merging
level = graythresh(img_diff);
bw = im2bw(img_diff,level) & mask_arena;
stats = regionprops(double(bw),'Centroid'); % get centroid of all white pixels
centroid = stats.Centroid;
imshowpair(img2,bw);


function edit_skip_frame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_skip_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_skip_frame as text
%        str2double(get(hObject,'String')) returns contents of edit_skip_frame as a double

handles.skip_frame = str2num(get(hObject,'String'));
display(sprintf('Now skipping every %d frames',handles.skip_frame));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_skip_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_skip_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ds = struct2cell2dataset(structure)
% searches for any char fields in a structure and converts them to cells
% then converts structure to dataset

fn = fieldnames(structure);
chars = structfun(@ischar, structure);
for n = find(chars)'
    structure.(fn{n}) = cellstr(structure.(fn{n}));
end

ds = struct2dataset(structure);


function edit_level_Callback(hObject, eventdata, handles)
% hObject    handle to edit_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_level as text
%        str2double(get(hObject,'String')) returns contents of edit_level as a double

handles.level = str2num(get(hObject,'String'));
display(sprintf('mouse threshold level set to %0.2f',handles.level));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_level_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_export.
function push_export_Callback(hObject, eventdata, handles)
% hObject    handle to push_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mds = getappdata(handles.figure1, 'mds');

if isempty(mds)
    beep
    display('no dataset found')
    return
end

[filename, pathname, filterspec] = uiputfile({'*.xlsx' 'EXCEL';'*.txt' 'Text'});

[~, fname, ext] = fileparts(filename);

if filterspec == 1
    run_text = false;
    
    try
        export(mds,'XLSfile',fullfile(pathname, [fname ext]))
    catch
        display('EXCEL export failed, saving as tab-delimited text')
        run_text = true;
    end
else
    run_text = true;
end

if run_text
    export(mds, 'file', fullfile(pathname, [fname '.txt']),'Delimiter', '\t')
end