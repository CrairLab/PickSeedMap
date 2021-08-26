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

% Last Modified by GUIDE v2.5 21-Feb-2020 17:39:41

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


if ~isempty(findobj('Tag', 'Manuvent_threshold'))
    MC_h = findobj('Tag', 'Manuvent_threshold');
    MC_data = guidata(MC_h);
    MC_passed = get(MC_data.Plot_correlation, 'UserData');
    plotCorrObj = MC_passed.plotCorrObj;

    %Get the current postition for the seed
    curPos = plotCorrObj.curPos;
    
    %Save the filename
    handles.Load_maps.UserData.filename = plotCorrObj.filename;

    clear MC_data MC_passed;

    %Plot the correlation map related to the seed
    corrM = CalculateCorrelationMap(plotCorrObj, handles);
    %handles.Status.Visible = 'On';  handles.Status.String = 'Finished!';

    %Calculate the averaged correlation near the seed
    if ~isfield(plotCorrObj, 'curTrace')
        avg_seed_corr = calculateSurrounding(curPos, corrM);
        handles.edit_seed.String = num2str(avg_seed_corr);
        handles.Save_data.UserData.avg_seed_corr = avg_seed_corr;
    end

    %Initiate Save_data.UserData
    handles.Save_data.UserData.save_strct = [];
    handles.Save_data.UserData.avg_rec_corr = [];
    handles.Save_data.UserData.hand_Position = [];
    %handles.Save_data.UserData.avg_hand_corr = [];
    %handles.Save_data.UserData.rec_Position = [];

    %Save the correlation matrix
    %handles.Load_maps.UserData.corrMatrix = corrM;
    handles.output.UserData.plotCorrObj = plotCorrObj;

    %Set default percentile value
    handles.Corr_percentile.UserData.curValue = 0.997;
    handles.Corr_percentile.String = num2str(0.997);
    disp('Default percentile for highly correlated region = 99.7%.')

    %Set the status text
    handles.Status.Visible = 'On';
    handles.Status.String = plotCorrObj.filename;
else
    disp('Did not load plotCorrObj from the Manuvent_threhold GUI!')
end


% UIWAIT makes PickSeedMap wait for user response (see UIRESUME)
% uiwait(handles.PickSeedMap);


function corrM = plotCorrelationMap(corrM, handles)
% Plot the correlation map w.r.t the current seed
    [y ,x] = findReccomandMax(corrM);    
    curPos = [x, y];    
    handles.Load_maps.UserData.curPos = curPos;
    
    %Specify current axis
    set(handles.PickSeedMap,'CurrentAxes',handles.CorrMap)
    
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

    
function corrM = CalculateCorrelationMap(plotCorrObj, handles)
% Plot the correlation map w.r.t the current seed
    A = plotCorrObj.curMovie;
    curPos = plotCorrObj.curPos;
    reg_flag = plotCorrObj.reg_flag;
    
    if reg_flag
        disp('Regress out the background!')
    else
        disp('Plot original correlation!')
    end
    
    %Get the fluorescent trace of the current roi
    if ~isfield(plotCorrObj, 'curTrace')
        seedTrace = A(curPos(2), curPos(1), :);
        seedTrace = seedTrace(:);
    else
        seedTrace = plotCorrObj.curTrace;
    end
    
    sz  = size(A);
    imgall = reshape(A, sz(1)*sz(2), sz(3));
    
    %Calculate correlation matrix
    try
        if ~reg_flag
            corrM = corr(imgall',seedTrace);
        else
            std_all = nanstd(imgall,0,2);
            lowPixels = std_all <= prctile(std_all,1);
            %avg_trace = nanmean(imgall,1);
            avg_trace = nanmean(imgall(lowPixels,:),1);
            corrM = partialcorr(imgall',seedTrace, avg_trace');
        end
        corrM = reshape(corrM, sz(1:2));

        %Plot the new correlation matrix
        corrM = plotCorrelationMap(corrM, handles);
    catch
        msgbox('Make sure traces have the same duration!')
    end

    
    
    
function [y2 ,x2] = findReccomandMax(corrM)
%Find the point that shows maximum correlation other than the seed given
%certain conditions

%Filter to mask out boundary artifact (due to rigid registration)
filter = ones(5);
corrM_conv = conv2(corrM,filter,'same');

%Find the maximum point
max_point = max(corrM_conv(:));

if max_point < 0
    %[y2 ,x2] = find(corrM_conv == max_point);
%else
    disp('Maximum correlation in the region is less then 0!')
    %warning('Try to find minimum point instead!')
    %min_point = min(corrM_conv(:));
    %[y2 ,x2] = find(corrM_conv == min_point);
end
[y2 ,x2] = find(corrM_conv == max_point);

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
    curFrame = handles.Frame.String;
    
    % Copy the obj to a new figure and save it
    f2 = figure;
    copyobj(handles.CorrMap, f2);
    set(f2.CurrentAxes, 'Units', 'Normalized');
    set(f2.CurrentAxes, 'OuterPosition', [0, 0, .95, .95]);
    colormap(jet);
    colorbar;
    
    saveas(f2, ['roi_', reg_flag, '_', curFrame, '_',...
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
    save(['roi_', num2str(reg_flag), '_', curFrame, '_',...
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
    
    %Update status
    handles.Status.Visible = 'On';
    handles.Status.String = file;
    
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
    
    %Clean old variables
    handles.output.UserData = [];
    handles.Save_data.UserData.save_strct = [];
    handles.Save_data.UserData.avg_rec_corr = [];
    handles.Save_data.UserData.hand_Position = [];

catch
    msgbox('Can not load correlation maps!','Error!')
end



function Frame_Callback(hObject, eventdata, handles)
% hObject    handle to Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Frame as text
%        str2double(get(hObject,'String')) returns contents of Frame as a double

try
    curIdx = str2double(get(hObject, 'String'));
    corrMatrix = handles.Load_maps.UserData.corrMatrix;
    corrM = corrMatrix(:,:,curIdx);
    plotCorrelationMap(corrM, handles);
    set(handles.Map_slider, 'Value', curIdx);
catch
    msgbox('Out of range. Can not jump to this map!', 'Error')
end

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

try
    corrMatrix = handles.Load_maps.UserData.corrMatrix;
    curIdx = round(get(hObject, 'Value'));
    corrM = corrMatrix(:,:,curIdx);
    plotCorrelationMap(corrM, handles);
    set(handles.Frame, 'String', num2str(curIdx));
catch
    msgbox('Wrong calling the slider callback, no corrMatrix loaded!', 'Error')
end

% --- Executes during object creation, after setting all properties.
function Map_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Map_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in Regionprops.
function Regionprops_Callback(hObject, eventdata, handles)
% hObject    handle to Regionprops (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    %corrM_masked = handles.Draw_region.UserData.corrM_masked;
    %Calculate the masked corrM
    BW = handles.Draw_region.UserData.BW;
    curMap = handles.CorrMap.UserData.curMap;     
    corrM_masked = curMap .* BW;
    handles.Draw_region.UserData.corrM_masked = corrM_masked;
    
    %Get the threshold
    threshold = str2double(handles.Regionprops_threshold.String);
    if isnan(threshold)
        warning('Input error! Please type in a number!')
        disp('Defulat threshold = 0.95!')
        threshold = 0.95;
    end
    
    %Area of the roi
    RegionArea = sum(BW(:));
    
    %Find the point that shows maximum correlation in a defined roi region
    HighCorrRegion = corrM_masked > threshold;
    figure; imshow(HighCorrRegion);
    stats = regionprops(HighCorrRegion, 'Area', 'Eccentricity', 'FilledArea',...
        'MajorAxisLength', 'MinorAxisLength', 'Orientation');
    AllArea = []; AllArea = [AllArea; stats.FilledArea];
    LargestRegion = AllArea == max(AllArea);
    LargestStats = stats(LargestRegion);
    
    %Normalized by the area of the roi
    if ~isempty(LargestStats)
        LargestStats.NormalizedArea = LargestStats.Area./RegionArea;
        LargestStats.NormalizedFilledArea = LargestStats.FilledArea./RegionArea;
        LargestStats.NormalizedMajorAxisLength = LargestStats.MajorAxisLength./RegionArea;
        LargestStats.NormalizedMinorAxisLength = LargestStats.MinorAxisLength./RegionArea;
        LargestStats.roiArea = RegionArea;
    else
        msgbox('Threshold too high. No pixel is above it!')
        return
    end
    
    %handles.Save_data.UserData.avg_rec_corr = avg_rec_corr;
    handles.Regionprops.UserData.LargestStats = LargestStats;
    handles.Status.Visible = 'On';
    handles.Status.String = 'Regionprops calculated!';
    
    %Save the regionprops result
    filename = handles.Load_maps.UserData.filename;
    %If the object was inherited from the Manuvent_threshold GUI
    if isfield(handles.output.UserData,'plotCorrObj')
        plotCorrObj = handles.output.UserData.plotCorrObj;
        reg_flag = plotCorrObj.reg_flag;
        disp('Loaded reg_flag!')
    else
        reg_flag = filename(20:21);
    end
    curPos = handles.Load_maps.UserData.curPos;
    uisave({'LargestStats'},['Regionprops_', num2str(reg_flag), '_', ...
        num2str(curPos(1)) '_',  num2str(curPos(2)),'.mat']);
    handles.Status.String = 'Regionprops saved!';
    disp(['Normalized filled area = ' num2str(LargestStats.NormalizedFilledArea)])
catch
    msgbox('Please define a roi region first!', 'Error');
end


function Regionprops_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to Regionprops_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Regionprops_threshold as text
%        str2double(get(hObject,'String')) returns contents of Regionprops_threshold as a double


% --- Executes during object creation, after setting all properties.
function Regionprops_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Regionprops_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Corr_region.
function Corr_region_Callback(hObject, eventdata, handles)
% hObject    handle to Corr_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try  
    %Get the current threshold
    try
        %Get and calculate the percentile
        curThreshold = str2double(get(handles.Corr_threshold, 'String'));
        handles.Corr_threshold.UserData.curValue = curThreshold;
        disp(['Plot region with correlation > ' num2str(curThreshold*100) '%'])
    catch
        msgbox('Please input a number btw 0-1!', 'Error!')
    end
    
    %Plot the highly correlated region
    curThreshold = handles.Corr_threshold.UserData.curValue;
    corrM = handles.CorrMap.UserData.curMap;
    %Renew correlation map first
    plotCorrelationMap(corrM, handles);
    [x,y] = find(corrM>=curThreshold);
    Correlated_region.x = x;
    Correlated_region.y = y;
    hObject.UserData.Correlated_region = Correlated_region;
    plot(handles.CorrMap, y,x,'w*')  
    
    %Quantify pixels with corr>threhold in each hemishere
    sz = size(corrM);
    mid = round(sz(2)/2);
    nLeft = sum(y<mid);
    nRight = sum(y>mid);
    ratio = nLeft/nRight;    
    diff = abs(nLeft - nRight);
    FOVsize = sum(~isnan(corrM(:)));
    savename = ['LR_balance' '_' num2str(curThreshold) '_' ...
        handles.Load_maps.UserData.filename];
    save(savename,'nLeft','nRight','ratio','diff','FOVsize')
catch
    msgbox('Can not get the correlated region!', 'Error!')
    return
end



function Corr_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to Corr_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Corr_threshold as text
%        str2double(get(hObject,'String')) returns contents of Corr_threshold as a double
try
    %Get and calculate the percentile
    curThreshold = str2double(get(hObject, 'String'));
    handles.Corr_threshold.UserData.curValue = curThreshold;
    disp(['Plot region with correlation > ' num2str(curThreshold*100) '%'])
catch
    msgbox('Please input a number btw 0-1!', 'Error!')
end



% --- Executes during object creation, after setting all properties.
function Corr_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Corr_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Plot_avg.
function Plot_avg_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_avg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


try
    %Get required data
    curThreshold = handles.Corr_threshold.UserData.curValue;
    corrM = handles.CorrMap.UserData.curMap;
    plotCorrObj = handles.output.UserData.plotCorrObj;
    curMovie = plotCorrObj.curMovie;
    
    %Reshape matrices
    sz = size(corrM);
    corrM = reshape(corrM, [sz(1)*sz(2),1]);
    curMovie = reshape(curMovie, [sz(1)*sz(2), size(curMovie,3)]);
    
    %Calculate the averaged trace
    Avg_trace = nanmean(curMovie(corrM>=curThreshold,:),1);
    hObject.UserData.Avg_trace = Avg_trace;
    
    %Plot the trace
    plot(handles.Trace, Avg_trace, 'LineWidth', 2);
    handles.Trace.XLim = [1, length(Avg_trace)];
    max_time = find(Avg_trace == max(Avg_trace));
    hold(handles.Trace, 'on');
    plot(handles.Trace, max_time, max(Avg_trace), 'r*')
    plot(handles.Trace, 10:14, Avg_trace(10:14), 'g', 'LineWidth', 2)
    Duration = length(Avg_trace);
    plot(handles.Trace, 1:Duration, 0.02*ones(Duration,1), 'r')
    
    %Save the values to UserData
    hObject.UserData.max_time = max_time;
    hObject.UserData.max_value = max(Avg_trace);
    
catch
    msgbox('Can not plot the averaged trace!', 'Error!')
    return
end

% --- Executes on button press in Save_trace.
function Save_trace_Callback(hObject, eventdata, handles)
% hObject    handle to Save_trace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curThreshold = handles.Corr_threshold.UserData.curValue;
Correlated_region = handles.Corr_region.UserData.Correlated_region;
Avg_trace = handles.Plot_avg.UserData.Avg_trace;
max_time = handles.Plot_avg.UserData.max_time;
max_value = handles.Plot_avg.UserData.max_value;
plotCorrObj = handles.output.UserData.plotCorrObj;
filename = plotCorrObj.filename;

%Locate specific tag
t = strfind(filename,'AveragedMatrix');
if ~isempty(t)
    filename = filename(t+15:end);
end

savename = [filename '_averaged_trace_' num2str(curThreshold)];
save([savename '.mat'], 'curThreshold', 'Correlated_region', 'Avg_trace', 'max_time', 'max_value')
saveas(handles.Trace, [savename '.png'])  


% --- Executes on button press in Clean_trace.
function Clean_trace_Callback(hObject, eventdata, handles)
% hObject    handle to Clean_trace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.Trace)


% --- Executes on button press in Regress_flag.
function Regress_flag_Callback(hObject, eventdata, handles)
% hObject    handle to Regress_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Regress_flag


% --- Executes on button press in Plot_correlation.
function Plot_correlation_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_correlation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    handles.Status.Visible = 'On';
    handles.Status.String = 'Select a pixel!';

    %Define roi
    set(handles.PickSeedMap,'CurrentAxes',handles.CorrMap)
    roi = drawpoint('Color', 'm');

    %Get roi position
    curPos = round(roi.Position); 

    %Update current roi position
    plotCorrObj = handles.output.UserData.plotCorrObj;
    plotCorrObj.curPos = curPos;
    
    %Update whether do background regression
    plotCorrObj.reg_flag = handles.Regress_flag.Value;
    
    %Plot the correlation map related to the seed
    CalculateCorrelationMap(plotCorrObj, handles);
catch
    msgbox('Something wrong, check if plotCorrObj correctly loaded!', 'Error')
end


% --- Executes on button press in ShowOldMask.
function ShowOldMask_Callback(hObject, eventdata, handles)
% hObject    handle to ShowOldMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
