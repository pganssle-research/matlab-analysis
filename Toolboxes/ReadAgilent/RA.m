function varargout = RA(varargin)

% RA M-file for RA.fig
%      RA, by itself, creates a new RA or raises the existing
%      singleton*.
%
%      H = RA returns the handle to a new RA or the handle to
%      the existing singleton*.
%
%      RA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RA.M with the given input arguments.
%
%      RA('Property','Value',...) creates a new RA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RV_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RA

% Last Modified by GUIDE v2.5 30-Mar-2011 13:25:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RA_OpeningFcn, ...
                   'gui_OutputFcn',  @RA_OutputFcn, ...
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


% --- Executes just before RA is made visible.
function RA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RA (see VARARGIN)

% Choose default command line output for RA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%fig = openfig(mfilename,'reuse');
   %a=imread('SimpsonGray.jpg');
   
 %  a=imread('VarianScreen.png');
  % b=imread('');
   %subplot(3,2,3)
   
  % H_image1=imshow(a);




% --------------------------------------------------------------------
function Varian_Open_Callback(hObject, eventdata, handles)
% hObject    handle to Varian_Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Open_kspace_Callback(hObject, eventdata, handles)
% hObject    handle to Open_kspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global image_space final_kspace
global pathname final_kspace
nslice=1;
[final_kspace,np,ntraces,nblocks]= ReadAgilentStart(nslice);

% --------------------------------------------------------------------
function Open_images_Callback(hObject, eventdata, handles)
% hObject    handle to Open_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear image_space dim rank1 bits1 
clear global image_space final_kspace nblocks array ntraces seqcon final_kspace
global pathname
ReadAgilentFDF



% --------------------------------------------------------------------
function Open_Spec_Callback(hObject, eventdata, handles)
% hObject    handle to Open_Spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Open_spec_Callback(hObject, eventdata, handles)
% hObject    handle to Open_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global final_kspace image_space np ntraces nblocks
global pathname final_kspace
nspec=1
[final_kspace,np,ntraces,nblocks]= ReadAgilentStartSpec(nspec);





