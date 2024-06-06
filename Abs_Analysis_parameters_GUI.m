function varargout = Abs_Analysis_parameters_GUI(varargin)
% Check if the window is already open
existingWindows = findobj('Type', 'figure', 'Name', 'Abs Analysis Parameters Loader', 'Visible', 'on');
if ~isempty(existingWindows)
    cprintf('Keywords','The Abs Analysis Parameters Loader window is already open.\n');
    return;
end

% Check if Abs_Analysis_parameters already exists in the workspace
if evalin('base', 'exist(''Abs_Analysis_parameters'', ''var'')')
    Abs_Analysis_parameters = evalin('base', 'Abs_Analysis_parameters');
else
    %% Initialize the GUI data
    %first and main ROI
    Abs_Analysis_parameters.roiColStart = 10;
    Abs_Analysis_parameters.roiColEnd = 2000;
    Abs_Analysis_parameters.roiRowStart = 10;
    Abs_Analysis_parameters.roiRowEnd = 2000;

    %second ROI
    Abs_Analysis_parameters.ROI2 = 0;
    Abs_Analysis_parameters.roi2ColStart = 650;
    Abs_Analysis_parameters.roi2ColEnd = 1450;
    Abs_Analysis_parameters.roi2RowStart= 1000;
    Abs_Analysis_parameters.roi2RowEnd = 1600;

    %third ROI
    Abs_Analysis_parameters.ROI3 = 0;
    Abs_Analysis_parameters.roi3ColStart = 800;
    Abs_Analysis_parameters.roi3ColEnd = 1200;
    Abs_Analysis_parameters.roi3RowStart = 700;
    Abs_Analysis_parameters.roi3RowEnd = 1500;

    %general parameters
    Abs_Analysis_parameters.roiStep = 10;
    Abs_Analysis_parameters.pixelmarkerX = 400;
    Abs_Analysis_parameters.pixelmarkerY = 800;
    Abs_Analysis_parameters.filter = 0;
    Abs_Analysis_parameters.camera = 'in-trap';
    Abs_Analysis_parameters.fittype = 'gauss2d';
    Abs_Analysis_parameters.MaxOD = 1;
    Abs_Analysis_parameters.roiZoom = 0;
    Abs_Analysis_parameters.Jointfit = 0;

    % Save data to the workspace
    assignin('base','Abs_Analysis_parameters',Abs_Analysis_parameters);
    cprintf('Keywords','Initialization of the imaging analysis parameters\n')
end

% Initialize the GUI and its components
h = figure('Visible','off','color','white','Position',[70,180,700,340], 'Resize', 'Off', 'Name', 'Abs Analysis Parameters Loader', 'NumberTitle', 'on', 'MenuBar', 'none', 'ToolBar', 'none');
set(h, 'WindowStyle', 'docked');

% Create labels and input fields
Font_Size = 15;
graycolor = [.32 .32 .32];
FrameColor = [0, 0, 0];
BackgroundColor2 = [0.99, 0.99, 0.99];
whitecolor = [1, 1, 1];
blackcolor = [0 0 0];
bluecolor = [0 0.339 0.839];
greencolor = [0.282 0.715 0.433];
goldcolor = [0.96 0.87 0.70];
redcolor = [.972 .081 .031];
Font_Name = 'Verdana';
adjustx = - 60;
adjusty = -370;
height = 25;
heigth_spacing = 7;
width_spacing = 10;
width_normal = 120;
width_double = 250;
%% Third ROI
uicontrol('Style','text','Position',[adjustx+70,adjusty+611+2*height+3*heigth_spacing-heigth_spacing,width_double,height],'String','roiRow 3','Tag', 'roiCol3Title','FontSize',Font_Size, 'BackgroundColor', blackcolor, 'ForegroundColor', whitecolor,'Visible', 'off', 'Tag', 'roiRow3Title','FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+70,adjusty+611+height+2*heigth_spacing-heigth_spacing,width_normal,height],'String',num2str(Abs_Analysis_parameters.roi3RowStart),'Tag','roi3RowStart','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'Visible', 'off','FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+200,adjusty+611+height+2*heigth_spacing-heigth_spacing,width_normal,height],'String',num2str(Abs_Analysis_parameters.roi3RowEnd),'Tag','roi3RowEnd','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'Visible', 'off','FontName', Font_Name);

uicontrol('Style','text','Position',[adjustx+330,adjusty+611+2*height+3*heigth_spacing-heigth_spacing,width_double,height],'String','roiCol 3','FontSize',Font_Size, 'BackgroundColor', blackcolor, 'ForegroundColor', whitecolor,'Visible', 'off', 'Tag', 'roiCol3Title','FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+330,adjusty+611+height+2*heigth_spacing-heigth_spacing,width_normal,height],'String',num2str(Abs_Analysis_parameters.roi3ColStart),'Tag','roi3ColStart','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'Visible', 'off','FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+460,adjusty+611+height+2*heigth_spacing-heigth_spacing,width_normal,height],'String',num2str(Abs_Analysis_parameters.roi3ColEnd),'Tag','roi3ColEnd','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'Visible', 'off','FontName', Font_Name);

% Create an activation toggle  button for third region of interest  with LED indicator
uicontrol('Style', 'togglebutton', 'String', 'Third ROI Off', 'Position', [adjustx+590,adjusty+574+height+42+heigth_spacing-heigth_spacing,240,height+42], ...
    'Tag', 'thirdRoiButton', 'Callback', {@ThirdROI_button_Callback,h},'FontSize',15, 'ForegroundColor', blackcolor, 'BackgroundColor', goldcolor,'Visible','off','FontName', Font_Name);

% Third Roi Toggle button callback function
    function ThirdROI_button_Callback(hObject, eventdata, h)
        % Check if the button is selected
        if (get(hObject, 'Value') == get(hObject, 'Max'))
            % Toggle button is pressed
            Abs_Analysis_parameters.ROI3 = 1;
            set(hObject, 'String', 'Third ROI On', 'FontSize', 15, 'BackgroundColor', bluecolor,'ForegroundColor',whitecolor,'FontName', Font_Name); % changing color to blue (on)
            set(findobj(h, 'Tag', 'roiRow3Title'), 'Visible', 'on');
            set(findobj(h, 'Tag', 'roi3RowStart'), 'Visible', 'on');
            set(findobj(h, 'Tag', 'roi3RowEnd'), 'Visible', 'on');
            set(findobj(h, 'Tag', 'roiCol3Title'), 'Visible', 'on');
            set(findobj(h, 'Tag', 'roi3ColStart'), 'Visible', 'on');
            set(findobj(h, 'Tag', 'roi3ColEnd'), 'Visible', 'on');
            set(findobj(h, 'Tag', 'thirdRoiButton'), 'Position', [adjustx+590,adjusty+574+height+42+heigth_spacing-heigth_spacing,240,height+42]);
        else
            % Toggle button is not pressed
            Abs_Analysis_parameters.ROI3 = 0;
            set(hObject, 'String', 'Third ROI Off', 'FontSize', 15, 'BackgroundColor', goldcolor,'ForegroundColor',blackcolor,'FontName', Font_Name); % changing color back to blue (off)
            set(findobj(h, 'Tag', 'roiRow3Title'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roi3RowStart'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roi3RowEnd'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roiCol3Title'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roi3ColStart'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roi3ColEnd'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'thirdRoiButton'), 'Position', [adjustx+70,adjusty+574+height+42+heigth_spacing-heigth_spacing,760,height+42]);
        end
    end

%% Second ROI
uicontrol('Style','text','Position',[adjustx+70,adjusty+540+2*height+3*heigth_spacing-heigth_spacing,width_double,height],'String','roiRow 2','Tag', 'roiCol2Title','FontSize',Font_Size, 'BackgroundColor', blackcolor, 'ForegroundColor', whitecolor,'Visible', 'off', 'Tag', 'roiRow2Title','FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+70,adjusty+540+height+2*heigth_spacing-heigth_spacing,width_normal,height],'String',num2str(Abs_Analysis_parameters.roi2RowStart),'Tag','roi2RowStart','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'Visible', 'off','FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+200,adjusty+540+height+2*heigth_spacing-heigth_spacing,width_normal,height],'String',num2str(Abs_Analysis_parameters.roi2RowEnd),'Tag','roi2RowEnd','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'Visible', 'off','FontName', Font_Name);

uicontrol('Style','text','Position',[adjustx+330,adjusty+540+2*height+3*heigth_spacing-heigth_spacing,width_double,height],'String','roiCol 2','FontSize',Font_Size, 'BackgroundColor', blackcolor, 'ForegroundColor', whitecolor,'Visible', 'off', 'Tag', 'roiCol2Title');
uicontrol('Style','edit','Position',[adjustx+330,adjusty+540+height+2*heigth_spacing-heigth_spacing,width_normal,height],'String',num2str(Abs_Analysis_parameters.roi2ColStart),'Tag','roi2ColStart','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'Visible', 'off','FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+460,adjusty+540+height+2*heigth_spacing-heigth_spacing,width_normal,height],'String',num2str(Abs_Analysis_parameters.roi2ColEnd),'Tag','roi2ColEnd','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'Visible', 'off','FontName', Font_Name);
% Create an activation toggle  button for second region of interest  with LED indicator
uicontrol('Style', 'togglebutton', 'String', 'Second ROI Off', 'Position', [adjustx+70,adjusty+500+height+42+heigth_spacing-heigth_spacing,760,height+42], ...
    'Tag', 'secondROIbutton', 'Callback', {@SecondROI_button_Callback,h},'FontSize',15, 'ForegroundColor', blackcolor, 'BackgroundColor',goldcolor,'FontName', Font_Name);
% Second Roi Toggle button callback function
    function SecondROI_button_Callback(hObject, eventdata, h)
        % Check if the button is selected
        if (get(hObject, 'Value') == get(hObject, 'Max'))
            % Toggle button is pressed
            Abs_Analysis_parameters.ROI2 = 1;
            set(hObject, 'String', 'Second ROI On', 'FontSize', 15, 'BackgroundColor', bluecolor,'ForegroundColor',whitecolor,'FontName', Font_Name); % changing color to blue (on)
            set(findobj(h, 'Tag', 'roiRow2Title'), 'Visible', 'on');
            set(findobj(h, 'Tag', 'roi2RowStart'), 'Visible', 'on');
            set(findobj(h, 'Tag', 'roi2RowEnd'), 'Visible', 'on');
            set(findobj(h, 'Tag', 'roiCol2Title'), 'Visible', 'on');
            set(findobj(h, 'Tag', 'roi2ColStart'), 'Visible', 'on');
            set(findobj(h, 'Tag', 'roi2ColEnd'), 'Visible', 'on');
            set(findobj(h, 'Tag', 'secondROIbutton'), 'Position', [adjustx+590,adjusty+500+height+42+heigth_spacing-heigth_spacing,240,height+42]);


            % Turn on the "third ROI" toggle button and make it dissapear

            thirdROIToggleButton = findobj('Tag', 'thirdRoiButton');
            set(thirdROIToggleButton, 'Value', get(thirdROIToggleButton, 'Min'));
            set(thirdROIToggleButton, 'String', 'Third ROI Off', 'FontSize', 15, 'BackgroundColor', goldcolor,'ForegroundColor',blackcolor,'FontName', Font_Name);
            set(findobj(h, 'Tag', 'thirdRoiButton'),'Position',[adjustx+70,adjusty+574+height+42+heigth_spacing-heigth_spacing,760,height+42],'Visible','On');
          
            %disable the Zoom buttong and set it ot off!
            zoomToggleButton = findobj(h, 'Tag', 'zoomButton');
             set(zoomToggleButton, 'Value', get(zoomToggleButton, 'Min'));
             set(zoomToggleButton, 'String', '<HTML>Z<BR>o<BR>o<BR>m<BR> <BR>O<BR>f<BR>f', 'FontSize', 12, 'BackgroundColor', [0.96 0.87 0.70],'FontName', Font_Name);
            set(zoomToggleButton, 'Enable', 'off');

        else
            % Toggle button is not pressed
            Abs_Analysis_parameters.ROI2 = 0;
            set(hObject, 'String', 'Second ROI Off', 'FontSize', 15, 'BackgroundColor', goldcolor,'ForegroundColor',blackcolor,'FontName', Font_Name); % changing color back to blue (off)
            set(findobj(h, 'Tag', 'roiRow2Title'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roi2RowStart'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roi2RowEnd'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roiCol2Title'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roi2ColStart'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roi2ColEnd'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'secondROIbutton'), 'Position', [adjustx+70,adjusty+500+height+42+heigth_spacing-heigth_spacing,760,height+42]);

            % Turn off the "third ROI" toggle button and make it dissapear

            thirdROIToggleButton = findobj('Tag', 'thirdRoiButton');
            set(thirdROIToggleButton, 'Value', get(thirdROIToggleButton, 'Min'));
            set(thirdROIToggleButton, 'String', 'Third ROI Off', 'FontSize', 15, 'BackgroundColor', [0.96 0.87 0.70],'FontName', Font_Name);
            set(findobj(h, 'Tag', 'thirdRoiButton'),'Visible','off');

            set(findobj(h, 'Tag', 'roiRow3Title'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roi3RowStart'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roi3RowEnd'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roiCol3Title'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roi3ColStart'), 'Visible', 'off');
            set(findobj(h, 'Tag', 'roi3ColEnd'), 'Visible', 'off');
            Abs_Analysis_parameters.ROI3 = 0;


            %Reenable the button zoom
            zoomToggleButton = findobj(h, 'Tag', 'zoomButton');
            set(zoomToggleButton, 'Enable', 'on');
        end
    end





%% Main ROI 
uicontrol('Style','text','Position',[adjustx+70,adjusty+540-heigth_spacing,width_double,height],'String','roiRow','FontSize',Font_Size, 'BackgroundColor', blackcolor, 'ForegroundColor', whitecolor,'FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+70,adjusty+500-heigth_spacing,width_normal,height],'String',num2str(Abs_Analysis_parameters.roiRowStart),'Tag','roiRowStart','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+200,adjusty+500-heigth_spacing,width_normal,height],'String',num2str(Abs_Analysis_parameters.roiRowEnd),'Tag','roiRowEnd','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'FontName', Font_Name);

uicontrol('Style','text','Position',[adjustx+330,adjusty+540-heigth_spacing,width_double,height],'String','roiCol','FontSize',Font_Size, 'BackgroundColor', blackcolor, 'ForegroundColor', whitecolor,'FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+330,adjusty+500-heigth_spacing,width_normal,height],'String',num2str(Abs_Analysis_parameters.roiColStart),'Tag','roiColStart','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+460,adjusty+500-heigth_spacing,width_normal,height],'String',num2str(Abs_Analysis_parameters.roiColEnd),'Tag','roiColEnd','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'FontName', Font_Name);

uicontrol('Style','text','Position',[adjustx+330,adjusty+460-heigth_spacing,width_normal,height],'String','roiStep','FontSize',Font_Size, 'BackgroundColor', blackcolor, 'ForegroundColor', whitecolor,'FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+330,adjusty+420,width_normal,height],'String',num2str(Abs_Analysis_parameters.roiStep),'Tag','roiStep','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'FontName', Font_Name);

uicontrol('Style','text','Position',[adjustx+70,adjusty+460-heigth_spacing,width_double,height],'String','Pixel Marker','FontSize',Font_Size, 'BackgroundColor', blackcolor, 'ForegroundColor', whitecolor,'FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+70,adjusty+420,width_normal,height],'String',num2str(Abs_Analysis_parameters.pixelmarkerX),'Tag','pixelmarkerX','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'FontName', Font_Name);
uicontrol('Style','edit','Position',[adjustx+200,adjusty+420,width_normal,height],'String',num2str(Abs_Analysis_parameters.pixelmarkerY),'Tag','pixelmarkerY','FontSize',Font_Size, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'FontName', Font_Name);

%% Create dropdown menus
%drop down Fittype
uicontrol('Style', 'text', 'Position', [adjustx+460,adjusty+460-heigth_spacing,width_normal,height], 'String', 'Fittype','FontSize',Font_Size, 'BackgroundColor', blackcolor, 'ForegroundColor', whitecolor,'FontName', Font_Name);
uicontrol('Style', 'popupmenu', 'String', 'gauss2d|twocomp2d', 'Position', [adjustx+460,adjusty+420+heigth_spacing/2,width_normal,height-3], 'Tag', 'fittype','FontSize',Font_Size-3, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'FontName', Font_Name);
%drop down Camera #
uicontrol('Style', 'text', 'Position', [adjustx+590,adjusty+460-heigth_spacing,120,height], 'String', 'Camera','FontSize',Font_Size, 'BackgroundColor', blackcolor, 'ForegroundColor', whitecolor,'FontName', Font_Name);
uicontrol('Style', 'popupmenu', 'String', 'in-trap|drop 2 ', 'Position', [adjustx+590,adjusty+420+heigth_spacing/2,width_normal,height-3], 'Tag', 'camera','FontSize',Font_Size-3, 'ForegroundColor', blackcolor, 'BackgroundColor', BackgroundColor2,'FontName', Font_Name);

%% Create a load last image button
uicontrol('Style', 'pushbutton', 'String', 'Load last image', 'Position', [218,5,197,35],'Callback', {@load_button_Callback,h},'FontSize',18, 'ForegroundColor', whitecolor, 'BackgroundColor', graycolor,'FontName', Font_Name);

%% Create a load image# button
uicontrol('Style', 'pushbutton', 'String', 'Load image#', 'Position', [428,5,168,35],'Callback', {@loadSpe_button_Callback,h},'FontSize',18, 'ForegroundColor', whitecolor, 'BackgroundColor', graycolor,'FontName', Font_Name);
%Create a input field for image number
uicontrol('Style', 'edit', 'String', '1', 'Position', [608,5,50,35], 'Tag', 'imageNum','FontSize',18, 'BackgroundColor', 'white', 'Callback', {@imageNum_input_Callback, h});

%% Create an update button
uicontrol('Style', 'pushbutton', 'String', 'Update', 'Position', [8,5,207,35],'Callback', {@update_button_Callback,h},'FontSize',18, 'ForegroundColor', whitecolor, 'BackgroundColor', graycolor,'FontName', Font_Name);

%% Create a filter toggle button with LED indicator
uicontrol('Style', 'togglebutton', 'String', 'Filter Off', 'Position', [adjustx+590,adjusty+500-heigth_spacing,120,height+42], 'Tag', 'filterButton', 'Callback', {@filter_button_Callback,h},'FontSize',14, 'ForegroundColor', whitecolor, 'BackgroundColor', redcolor,'FontName', Font_Name);

%% Create a Full frame button button
fullFrameButton = uicontrol('Style', 'pushbutton', 'String', '<HTML>F<BR>u<BR>l<BR>l<BR>f<BR>r<BR>a<BR>m<BR>e',...
    'Position', [adjustx+764,adjusty+380-heigth_spacing+5,30,height+158],...
    'Callback', {@fullframe_button_Callback,h} ,'FontSize',10, 'ForegroundColor', whitecolor, 'BackgroundColor', graycolor);

%% Create a Zoom toggle button with LED indicator
zoomButton = uicontrol('Style', 'togglebutton', 'String', '<HTML>Z<BR>o<BR>o<BR>m<BR> <BR>O<BR>f<BR>f', 'FontSize',12,...
    'Position', [adjustx+804,adjusty+404-heigth_spacing-18,30,height+158], 'Tag', 'zoomButton', 'Callback', {@Zoom_button_Callback,h},'FontSize',11, 'ForegroundColor', blackcolor, 'BackgroundColor',  [0.96 0.87 0.70],'FontName', Font_Name); % initial color set to blue (off)


%% Create a vertical slider for maxOD
uicontrol('Style','text','Position',[adjustx+720,adjusty+535,35,25],'String','OD','FontSize',Font_Size-4, 'BackgroundColor', blackcolor, 'ForegroundColor', whitecolor,'FontName', Font_Name);

minMaxOD = [0.05, 5]; % Minimum and maximum values for maxOD
stepSize = [0.1, 0.25]; % Step sizes for the slider
sliderPosition = [adjustx+725,adjusty+405,30,125];
sliderValue = min(max(Abs_Analysis_parameters.MaxOD, minMaxOD(1)), minMaxOD(2));

slider_axes = axes('Position', [sliderPosition(1)/h.Position(3), sliderPosition(2)/h.Position(4), sliderPosition(3)/h.Position(3), sliderPosition(4)/h.Position(4)]);
slider_axes.Visible = 'off';

sld = uicontrol('Style', 'slider', 'Parent', h, 'Units', 'pixels', 'Position', sliderPosition, 'Value', sliderValue, 'Min', minMaxOD(1), 'Max', minMaxOD(2), 'SliderStep', stepSize./(minMaxOD(2)-minMaxOD(1)), 'Callback', @slider_callback, 'Tag','ODSlider');

slider_label = uicontrol('Style', 'edit', 'Parent', h, 'Units', 'pixels', 'Position', [sliderPosition(1), sliderPosition(2)-30, sliderPosition(3), 20], 'String', num2str(sliderValue), 'HorizontalAlignment', 'right', 'FontSize', 12, 'Callback', @edit_callback);
%__________________________Callback functions
%% OD slider Callback
% Slider callback function
    function slider_callback(src, ~)
        value = src.Value;
        slider_label.String = num2str(value);
        h = src.Parent; % get the handle of the figure
        Abs_Analysis_parameters = getappdata(h, 'Abs_Analysis_parameters'); % load the data
        Abs_Analysis_parameters.MaxOD = value; % modify the value
        setappdata(h, 'Abs_Analysis_parameters', Abs_Analysis_parameters); % save the modified data
    end

% Edit OD field callback function
    function edit_callback(src, ~)
        value = str2double(src.String);
        if isnan(value) || value < minMaxOD(1) || value > minMaxOD(2)
            value = sld.Value; % If input is invalid, revert to slider value
        end
        src.String = num2str(value); % Ensure string is properly formatted
        sld.Value = value; % Update slider value
        h = src.Parent; % get the handle of the figure
        Abs_Analysis_parameters = getappdata(h, 'Abs_Analysis_parameters'); % load the data
        Abs_Analysis_parameters.MaxOD = value; % modify the value
        setappdata(h, 'Abs_Analysis_parameters', Abs_Analysis_parameters); % save the modified data
    end



% Make figure visible after adding all components
set(h,'Visible','on');

%% Load button callback function
    function load_button_Callback(hObject, eventdata, h)


        Abs_Analysis_parameters.roiColStart = str2double(get(findobj(h, 'Tag', 'roiColStart'),'String'));
        Abs_Analysis_parameters.roiColEnd = str2double(get(findobj(h, 'Tag', 'roiColEnd'),'String'));
        Abs_Analysis_parameters.roiRowStart = str2double(get(findobj(h, 'Tag', 'roiRowStart'),'String'));
        Abs_Analysis_parameters.roiRowEnd = str2double(get(findobj(h, 'Tag', 'roiRowEnd'),'String'));
        Abs_Analysis_parameters.roiStep = str2double(get(findobj(h, 'Tag', 'roiStep'),'String'));
        Abs_Analysis_parameters.pixelmarkerX = str2double(get(findobj(h, 'Tag', 'pixelmarkerX'),'String'));
        Abs_Analysis_parameters.pixelmarkerY = str2double(get(findobj(h, 'Tag', 'pixelmarkerY'),'String'));
        Abs_Analysis_parameters.filter = get(findobj(h, 'Tag', 'filterButton'),'Value');
        fittype_options = {'gauss2d','twocomp2d'};
        Abs_Analysis_parameters.fittype = fittype_options{get(findobj(h, 'Tag', 'fittype'),'Value')};
        camera_options = {'in-trap', 'drop 2'};
        Abs_Analysis_parameters.camera = camera_options{get(findobj(h, 'Tag', 'camera'),'Value')};
        %     Abs_Analysis_parameters.MaxOD = getappdata(h, 'Abs_Analysis_parameters').MaxOD;
        Abs_Analysis_parameters.MaxOD = get(findobj(h, 'Tag', 'ODSlider'),'Value');

        Abs_Analysis_parameters.roiZoom = get(findobj(h, 'Tag', 'zoomButton'),'Value');

                    Abs_Analysis_parameters.Jointfit = 0;%need to fix this soon

        Abs_Analysis_parameters.roi2ColStart = str2double(get(findobj(h, 'Tag', 'roi2ColStart'),'String'));
        Abs_Analysis_parameters.roi2ColEnd = str2double(get(findobj(h, 'Tag', 'roi2ColEnd'),'String'));
        Abs_Analysis_parameters.roi2RowStart = str2double(get(findobj(h, 'Tag', 'roi2RowStart'),'String'));
        Abs_Analysis_parameters.roi2RowEnd = str2double(get(findobj(h, 'Tag', 'roi2RowEnd'),'String'));
        Abs_Analysis_parameters.ROI2 = get(findobj(h, 'Tag', 'secondROIbutton'),'Value');

        Abs_Analysis_parameters.roi3ColStart = str2double(get(findobj(h, 'Tag', 'roi3ColStart'),'String'));
        Abs_Analysis_parameters.roi3ColEnd = str2double(get(findobj(h, 'Tag', 'roi3ColEnd'),'String'));
        Abs_Analysis_parameters.roi3RowStart = str2double(get(findobj(h, 'Tag', 'roi3RowStart'),'String'));
        Abs_Analysis_parameters.roi3RowEnd = str2double(get(findobj(h, 'Tag', 'roi3RowEnd'),'String'));
        Abs_Analysis_parameters.ROI3 = get(findobj(h, 'Tag', 'secondROIbutton'),'Value');
        
        % Save data to the workspace
        assignin('base','Abs_Analysis_parameters',Abs_Analysis_parameters);
        cprintf('Keywords','Loading last image with updated parameters\n');
        Abs_Analysis_GUI();
    end


%% Load a specific image button callback function
    function loadSpe_button_Callback(hObject, eventdata, h)


        Abs_Analysis_parameters.roiColStart = str2double(get(findobj(h, 'Tag', 'roiColStart'),'String'));
        Abs_Analysis_parameters.roiColEnd = str2double(get(findobj(h, 'Tag', 'roiColEnd'),'String'));
        Abs_Analysis_parameters.roiRowStart = str2double(get(findobj(h, 'Tag', 'roiRowStart'),'String'));
        Abs_Analysis_parameters.roiRowEnd = str2double(get(findobj(h, 'Tag', 'roiRowEnd'),'String'));
        Abs_Analysis_parameters.roiStep = str2double(get(findobj(h, 'Tag', 'roiStep'),'String'));
        Abs_Analysis_parameters.pixelmarkerX = str2double(get(findobj(h, 'Tag', 'pixelmarkerX'),'String'));
        Abs_Analysis_parameters.pixelmarkerY = str2double(get(findobj(h, 'Tag', 'pixelmarkerY'),'String'));
        Abs_Analysis_parameters.filter = get(findobj(h, 'Tag', 'filterButton'),'Value');
        fittype_options = {'gauss2d','twocomp2d'};
        Abs_Analysis_parameters.fittype = fittype_options{get(findobj(h, 'Tag', 'fittype'),'Value')};
        camera_options = {'in-trap', 'drop 2'};
        Abs_Analysis_parameters.camera = camera_options{get(findobj(h, 'Tag', 'camera'),'Value')};
        Abs_Analysis_parameters.MaxOD = get(findobj(h, 'Tag', 'ODSlider'),'Value');
        Abs_Analysis_parameters.roiZoom = get(findobj(h, 'Tag', 'zoomButton'),'Value');

            Abs_Analysis_parameters.Jointfit = 0;%need to fix this soon


        Abs_Analysis_parameters.roi2ColStart = str2double(get(findobj(h, 'Tag', 'roi2ColStart'),'String'));
        Abs_Analysis_parameters.roi2ColEnd = str2double(get(findobj(h, 'Tag', 'roi2ColEnd'),'String'));
        Abs_Analysis_parameters.roi2RowStart = str2double(get(findobj(h, 'Tag', 'roi2RowStart'),'String'));
        Abs_Analysis_parameters.roi2RowEnd = str2double(get(findobj(h, 'Tag', 'roi2RowEnd'),'String'));
        Abs_Analysis_parameters.ROI2 = get(findobj(h, 'Tag', 'secondROIbutton'),'Value');

        Abs_Analysis_parameters.roi3ColStart = str2double(get(findobj(h, 'Tag', 'roi3ColStart'),'String'));
        Abs_Analysis_parameters.roi3ColEnd = str2double(get(findobj(h, 'Tag', 'roi3ColEnd'),'String'));
        Abs_Analysis_parameters.roi3RowStart = str2double(get(findobj(h, 'Tag', 'roi3RowStart'),'String'));
        Abs_Analysis_parameters.roi3RowEnd = str2double(get(findobj(h, 'Tag', 'roi3RowEnd'),'String'));
        Abs_Analysis_parameters.ROI3 = get(findobj(h, 'Tag', 'secondROIbutton'),'Value');

        imageNum = str2double(get(findobj(h, 'Tag', 'imageNum'),'String'));
        
        if isnan(imageNum) || imageNum < 1
            cprintf('Errors','Invalid image number.\n');
            return;
        end
        imageNum = floor(imageNum); % ensure integer

        % Save data to the workspace
        assignin('base','Abs_Analysis_parameters',Abs_Analysis_parameters);
        cprintf('Keywords', sprintf('Loading image #%d\n', imageNum));
        Abs_Analysis_GUI('last', imageNum);
    end



% imageNum input field callback function
    function imageNum_input_Callback(hObject, eventdata, h)
        imageNum = str2double(get(hObject,'String'));
        if isnan(imageNum) || imageNum < 1
            set(hObject, 'String', '1'); % reset to default
            cprintf('Errors','Invalid image number. Please input a positive integer.\n');
        else
            set(hObject, 'String', num2str(floor(imageNum))); % ensure integer
        end
    end



%% update parameter specific image button callback function
    function update_button_Callback(hObject, eventdata, h)

        Abs_Analysis_parameters.roiColStart = str2double(get(findobj(h, 'Tag', 'roiColStart'),'String'));
        Abs_Analysis_parameters.roiColEnd = str2double(get(findobj(h, 'Tag', 'roiColEnd'),'String'));
        Abs_Analysis_parameters.roiRowStart = str2double(get(findobj(h, 'Tag', 'roiRowStart'),'String'));
        Abs_Analysis_parameters.roiRowEnd = str2double(get(findobj(h, 'Tag', 'roiRowEnd'),'String'));
        Abs_Analysis_parameters.roiStep = str2double(get(findobj(h, 'Tag', 'roiStep'),'String'));
        Abs_Analysis_parameters.pixelmarkerX = str2double(get(findobj(h, 'Tag', 'pixelmarkerX'),'String'));
        Abs_Analysis_parameters.pixelmarkerY = str2double(get(findobj(h, 'Tag', 'pixelmarkerY'),'String'));
        Abs_Analysis_parameters.filter = get(findobj(h, 'Tag', 'filterButton'),'Value');
        fittype_options = {'gauss2d','twocomp2d'};
        Abs_Analysis_parameters.fittype = fittype_options{get(findobj(h, 'Tag', 'fittype'),'Value')};
        camera_options = {'in-trap', 'drop 2'};
        Abs_Analysis_parameters.camera = camera_options{get(findobj(h, 'Tag', 'camera'),'Value')};
        Abs_Analysis_parameters.MaxOD = get(findobj(h, 'Tag', 'ODSlider'),'Value');
        Abs_Analysis_parameters.roiZoom = get(findobj(h, 'Tag', 'zoomButton'),'Value');

                    Abs_Analysis_parameters.Jointfit = 0;%need to fix this soon

        Abs_Analysis_parameters.roi2ColStart = str2double(get(findobj(h, 'Tag', 'roi2ColStart'),'String'));
        Abs_Analysis_parameters.roi2ColEnd = str2double(get(findobj(h, 'Tag', 'roi2ColEnd'),'String'));
        Abs_Analysis_parameters.roi2RowStart = str2double(get(findobj(h, 'Tag', 'roi2RowStart'),'String'));
        Abs_Analysis_parameters.roi2RowEnd = str2double(get(findobj(h, 'Tag', 'roi2RowEnd'),'String'));
        Abs_Analysis_parameters.ROI2 = get(findobj(h, 'Tag', 'secondROIbutton'),'Value');

          Abs_Analysis_parameters.roi3ColStart = str2double(get(findobj(h, 'Tag', 'roi3ColStart'),'String'));
        Abs_Analysis_parameters.roi3ColEnd = str2double(get(findobj(h, 'Tag', 'roi3ColEnd'),'String'));
        Abs_Analysis_parameters.roi3RowStart = str2double(get(findobj(h, 'Tag', 'roi3RowStart'),'String'));
        Abs_Analysis_parameters.roi3RowEnd = str2double(get(findobj(h, 'Tag', 'roi3RowEnd'),'String'));
        Abs_Analysis_parameters.ROI3 = get(findobj(h, 'Tag', 'secondROIbutton'),'Value');

        % Save data to the workspace
        assignin('base','Abs_Analysis_parameters',Abs_Analysis_parameters);
        cprintf('Keywords','New image analysis parameters uploaded\n');
    end


%% Full Frame button callback function
    function fullframe_button_Callback(hObject, eventdata, h)
        Abs_Analysis_parameters.roiRowStart = 10;
        Abs_Analysis_parameters.roiRowEnd = 2000;
        Abs_Analysis_parameters.roiColStart = 10;
        Abs_Analysis_parameters.roiColEnd = 2000;

        % Update the UI fields
        set(findobj(h, 'Tag', 'roiRowStart'), 'String', num2str(Abs_Analysis_parameters.roiRowStart));
        set(findobj(h, 'Tag', 'roiRowEnd'), 'String', num2str(Abs_Analysis_parameters.roiRowEnd));
        set(findobj(h, 'Tag', 'roiColStart'), 'String', num2str(Abs_Analysis_parameters.roiColStart));
        set(findobj(h, 'Tag', 'roiColEnd'), 'String', num2str(Abs_Analysis_parameters.roiColEnd));

        % Turn off the "Second ROI" toggle button
        secondROIToggleButton = findobj('Tag', 'secondROIbutton');
        set(secondROIToggleButton, 'Value', get(secondROIToggleButton, 'Min'));
        set(secondROIToggleButton, 'String', 'Second ROI Off', 'FontSize', 15, 'BackgroundColor', goldcolor,'ForegroundColor',blackcolor,'FontName', Font_Name);
        set(findobj(h, 'Tag', 'secondROIbutton'), 'Position', [adjustx+70,adjusty+500+height+42+heigth_spacing-heigth_spacing,760,height+42]);

        set(findobj(h, 'Tag', 'roiRow2Title'), 'Visible', 'off');
        set(findobj(h, 'Tag', 'roi2RowStart'), 'Visible', 'off');
        set(findobj(h, 'Tag', 'roi2RowEnd'), 'Visible', 'off');
        set(findobj(h, 'Tag', 'roiCol2Title'), 'Visible', 'off');
        set(findobj(h, 'Tag', 'roi2ColStart'), 'Visible', 'off');
        set(findobj(h, 'Tag', 'roi2ColEnd'), 'Visible', 'off');
        Abs_Analysis_parameters.ROI2 = 0;

        % Turn off the "third ROI" toggle button and make it dissapear

        thirdROIToggleButton = findobj('Tag', 'thirdRoiButton');
        set(thirdROIToggleButton, 'Value', get(thirdROIToggleButton, 'Min'));
        set(thirdROIToggleButton, 'String', 'Third ROI Off', 'FontSize', 15, 'BackgroundColor', [0.96 0.87 0.70],'FontName', Font_Name);
        set(findobj(h, 'Tag', 'thirdRoiButton'),'Visible','off');

        set(findobj(h, 'Tag', 'roiRow3Title'), 'Visible', 'off');
        set(findobj(h, 'Tag', 'roi3RowStart'), 'Visible', 'off');
        set(findobj(h, 'Tag', 'roi3RowEnd'), 'Visible', 'off');
        set(findobj(h, 'Tag', 'roiCol3Title'), 'Visible', 'off');
        set(findobj(h, 'Tag', 'roi3ColStart'), 'Visible', 'off');
        set(findobj(h, 'Tag', 'roi3ColEnd'), 'Visible', 'off');
        Abs_Analysis_parameters.ROI3 = 0;

        %Reenable the Zoom button
        zoomToggleButton = findobj(h, 'Tag', 'zoomButton');
        set(zoomToggleButton, 'Enable', 'on');


        setappdata(h, 'Abs_Analysis_parameters', Abs_Analysis_parameters);


    end


% Filter button callback function
    function filter_button_Callback(hObject, eventdata, h)
        state = get(hObject,'Value');
        if state
            set(hObject, 'String', 'Filter On', 'BackgroundColor', greencolor,'FontName', Font_Name);
        else
            set(hObject, 'String', 'Filter Off', 'BackgroundColor', redcolor,'FontName', Font_Name);
        end
    end



% Zoom button callback function
function Zoom_button_Callback(hObject, eventdata, h)
% Check if the button is selected
if (get(hObject,'Value') == get(hObject,'Max'))
    % Toggle button is pressed
    Abs_Analysis_parameters.roiZoom = 1;
    set(hObject, 'String', '<HTML>Z<BR>o<BR>o<BR>m<BR> <BR>O<BR>n','FontSize',11, 'BackgroundColor', bluecolor,'ForegroundColor',whitecolor,'FontName', Font_Name); % changing color to blue (on)
else
    % Toggle button is not pressed
    Abs_Analysis_parameters.roiZoom = 0;
    set(hObject, 'String', '<HTML>Z<BR>o<BR>o<BR>m<BR> <BR>O<BR>f<BR>f<BR>', 'FontSize',11,'BackgroundColor', goldcolor,'ForegroundColor',blackcolor,'FontName', Font_Name); % changing color back to blue (off)
end
end
end
