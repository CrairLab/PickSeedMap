function varargout = PickSeedMap(varargin)
% PICKSEEDMAP MATLAB code for PickSeedMap.fig
%      PICKSEEDMAP, by itself, creates a new PICKSEEDMAP or raises the existing
%      singleton*.
%
%      H = PICKSEEDMAP returns the handle to a new PICKSEEDMAP or the handle to
%      the existing singleton*.
%
%      PICKSEEDMAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PICKSEEDMAP.M with the given input arguments.
%
%      PICKSEEDMAP('Property','Value',...) creates a new PICKSEEDMAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PickSeedMap_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PickSeedMap_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PickSeedMap

% Last Modified by GUIDE v2.5 26-Dec-2019 15:14:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PickSeedMap_OpeningFcn, ...
                   'gui_OutputFcn',  @PickSeedMap_OutputFcn, ...
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


% --- Executes just before PickSeedMap is made visible.
function PickSeedMap_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PickSeedMap (see VARARGIN)

% Choose default command line output for PickSeedMap
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

%Initiate Save_data.UserData
%handles.Save_data.UserData.avg_rec_corr = [];
handles.Save_data.UserData.avg_hand_corr = [];
handles.Save_data.UserData.save_strct = [];
handles.Save_data.UserData.hand_Position = [];


% UIWAIT makes PickSeedMap wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function corrM = plotCorrelationMap(corrM, handles)
% Plot the correlation map w.r.t the current seed
    [y ,x] = findReccomandMax(corrM);    
    curPos = [x, y];    
    handles.Load_maps.UserData.curPos = curPos;
    
    %Plot the correlation map
    hold off;
    im = imagesc(handles.CorrMap, corrM); colormap jet; colorbar; axis image
    caxis([-0.2, 1]); title(['roi:', num2str(curPos)]);
    
    %Save current map
    handles.CorrMap.UserData.curMap = corrM;
   
    %Allow roi selection
    set(im, 'ButtonDownFcn', {@markEvents, handles});
    hold on
    
    %Label the seed
    x1 = curPos(1);
    y1 = curPos(2);
    fill([x1-2,x1-2,x1+2,x1+2],[y1-2,y1+2,y1+2,y1-2], 'y')
    
    %Calculate the averaged correlation near the seed
    avg_seed_corr = calculateSurrounding(curPos, corrM);
    handles.edit_seed.String = num2str(avg_seed_corr);
    handles.Save_data.UserData.avg_seed_corr = avg_seed_corr;

    
function [y2 ,x2] = findReccomandMax(corrM)
%Find the point that shows maximum correlation other than the seed given
%certain conditions

%Filter to mask out boundary artifact (due to rigid registration)
filter = ones(5);
corrM_conv = conv2(corrM,filter,'same');

%Find the maximum point
max_point = max(corrM_conv(:));

if max_point > 0
    [y2 ,x2] = find(corrM_conv == max_point);
else
    warning('Maximum correlation in the region is less then 0!')
    warning('Try to find minimum point instead!')
    min_point = min(corrM_conv(:));
    [y2 ,x2] = find(corrM_conv == min_point);
end


function save_strct = Get_stat(corrM)
%Get statistcs of the input correlation matrix
corrM(corrM == 0) = nan;

%Get statistics
avg_corr = nanmean(corrM(:)); %get the averaged correlation
median_corr = nanmedian(corrM(:)); %get the median correlation
[y2 ,x2] = findReccomandMax(corrM); %find the point with max correlation
max_corr = corrM(y2,x2); %Get the max correlation

%Save statistics
save_strct.avg_corr = avg_corr;
save_strct.median_corr = median_corr;
save_strct.max_position = [x2, y2];
save_strct.max_corr = max_corr;
save_strct.corrM = corrM;

try
    fill([x2-2,x2-2,x2+2,x2+2],[y2-2,y2+2,y2+2,y2-2], 'm')
catch
    msgbox('Detect more than one maximums/minimums!','Error')
end

function avg_corr = calculateSurrounding(curPos, corrM)
    x = curPos(1);
    y = curPos(2);
    surrond = corrM(y-2:y+2,x-2:x+2);
    avg_corr = nanmean(surrond(:));

function markEvents(h,~,handles)
%This function will allow user to mark a new event as well as store
%information of the defined roi (including the postition, the frame indices
%when the roi/event was created/initiated and deleted/ended)
%h        handle of the current image
%handles  handles of the GUI
roi = drawpoint(h.Parent, 'Color', 'g'); %Drawing a new roi
curPos = round(roi.Position);  %Current xy coordinates
corrM = handles.CorrMap.UserData.curMap;
avg_hand_corr = calculateSurrounding(curPos, corrM);
handles.edit_hand.String = num2str(avg_hand_corr);
handles.Save_data.UserData.avg_hand_corr = avg_hand_corr;
handles.Save_data.UserData.hand_Position = curPos;

%Listening to the moving events
addlistener(roi, 'ROIMoved', @(src,evt)movedCallback(src,evt,handles));


function movedCallback(roi,~,handles)
%Actions to take when moved the roi
%roi     the Point obj
%handles     handles of the GUI
curPos = round(roi.Position);  %Current xy coordinates
corrM = handles.CorrMap.UserData.curMap;
avg_hand_corr = calculateSurrounding(curPos, corrM);
handles.edit_hand.String = num2str(avg_hand_corr);
handles.Save_data.UserData.avg_hand_corr = avg_hand_corr;
handles.Save_data.UserData.hand_Position = curPos;


% --- Outputs from this function are returned to the command line.
function varargout = PickSeedMap_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_seed_Callback(hObject, eventdata, handles)
% hObject    handle to edit_seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_seed as text
%        str2double(get(hObject,'String')) returns contents of edit_seed as a double


% --- Executes during object creation, after setting all properties.
function edit_seed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_rec_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rec as text
%        str2double(get(hObject,'String')) returns contents of edit_rec as a double


% --- Executes during object creation, after setting all properties.
function edit_rec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_hand_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hand as text
%        str2double(get(hObject,'String')) returns contents of edit_hand as a double


% --- Executes during object creation, after setting all properties.
function edit_hand_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Save_data.
function Save_data_Callback(hObject, eventdata, handles)
% hObject    handle to Save_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    filename = handles.Load_maps.UserData.filename;
    reg_flag = filename(20:21);
    curPos = handles.Load_maps.UserData.curPos;
    saveas(handles.CorrMap, ['roi_', reg_flag, '_',...
        num2str(curPos(1)) '_',  num2str(curPos(2)),'.png'])  
    %avg_seed_corr = handles.Save_data.UserData.avg_seed_corr;
    %avg_rec_corr = handles.Save_data.UserData.avg_rec_corr;
    %avg_hand_corr = handles.Save_data.UserData.avg_hand_corr;
    %hand_Position = handles.Save_data.UserData.hand_Position;
    save_strct = handles.Save_data.UserData.save_strct;
    %save(['roi_', num2str(reg_flag), '_',...
        %num2str(curPos(1)) '_',  num2str(curPos(2)),'.mat'],...
        %'avg_seed_corr', 'avg_rec_corr', 'avg_hand_corr', 'curPos',...
        %'save_strct', 'hand_Position');       
    save(['roi_', num2str(reg_flag), '_', ...
        num2str(curPos(1)) '_',  num2str(curPos(2)),'.mat'], 'save_strct');
catch
    warning('Variables not defined!')
end


% --- Executes on button press in Draw_region.
function Draw_region_Callback(hObject, eventdata, handles)
% hObject    handle to Draw_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Status.Visible = 'On';
handles.Status.String = 'Pls draw a region!';
BW = roipoly;
handles.Draw_region.UserData.BW = BW;
handles.Save_data.UserData.save_strct.BW = BW;
corrM = handles.CorrMap.UserData.curMap;
corrM_masked = corrM .* BW;
handles.Draw_region.UserData.corrM_masked = corrM_masked;
handles.Status.String = 'ROI defined!';


% --- Executes on button press in Get_statistics.
function Get_statistics_Callback(hObject, eventdata, handles)
% hObject    handle to Get_statistics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    %corrM_masked = handles.Draw_region.UserData.corrM_masked;
    %Calculate the masked corrM
    BW = handles.Draw_region.UserData.BW;
    curMap = handles.CorrMap.UserData.curMap;     
    corrM_masked = curMap .* BW;
    handles.Draw_region.UserData.corrM_masked = corrM_masked;
    
    %Find the point that shows maximum correlation in a defined roi region
    save_strct = Get_stat(corrM_masked);
    %avg_rec_corr = save_strct.avg_corr;
    max_corr = save_strct.max_corr;
    handles.edit_rec.String = num2str(max_corr);
    %handles.Save_data.UserData.avg_rec_corr = avg_rec_corr;
    handles.Save_data.UserData.save_strct = save_strct;
    handles.Status.Visible = 'On';
    handles.Status.String = 'Statistics calculated!';
catch
    msgbox('Please define a roi region first!', 'Error');
end


% --- Executes during object creation, after setting all properties.
function Save_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Save_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in Load_maps.
function Load_maps_Callback(hObject, eventdata, handles)
% hObject    handle to Load_maps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Load the dF/F .mat movie
try
    [file,path]=uigetfile('*.mat','Please load the correlation matrix (mat) file!');
    load(fullfile(path,file));
    hObject.UserData.corrMatrix = corrMatrix;
    %hObject.UserData.filename = file;
    hObject.UserData.filename = file;
    sz = size(corrMatrix);
    
    %Set edit text
    handles.Frame.String = '1';
    
    %Set parameters for Map_slider
    set(handles.Map_slider, 'Min', 1);
    set(handles.Map_slider, 'Max', sz(3));
    set(handles.Map_slider, 'Value', 1);
    set(handles.Map_slider, 'SliderStep', [1/(sz(3)-1), 0.05]);
    
    %Plot the first map
    corrM = corrMatrix(:,:,1);
    plotCorrelationMap(corrM, handles);

catch
    msgbox('Can not load correlation maps!','Error!')
end



function Frame_Callback(hObject, eventdata, handles)
% hObject    handle to Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Frame as text
%        str2double(get(hObject,'String')) returns contents of Frame as a double

curIdx = str2double(get(hObject, 'String'));
corrMatrix = handles.Load_maps.UserData.corrMatrix;
corrM = corrMatrix(:,:,curIdx);
plotCorrelationMap(corrM, handles);
set(handles.Map_slider, 'Value', curIdx);


% --- Executes during object creation, after setting all properties.
function Frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function Map_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Map_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

corrMatrix = handles.Load_maps.UserData.corrMatrix;
curIdx = round(get(hObject, 'Value'));
corrM = corrMatrix(:,:,curIdx);
plotCorrelationMap(corrM, handles);
set(handles.Frame, 'String', num2str(curIdx));


% --- Executes during object creation, after setting all properties.
function Map_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Map_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
