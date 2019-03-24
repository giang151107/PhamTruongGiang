function varargout = Audio_Sign_Processing_in_matlab(varargin)
% AUDIO_SIGN_PROCESSING_IN_MATLAB MATLAB code for Audio_Sign_Processing_in_matlab.fig
%      AUDIO_SIGN_PROCESSING_IN_MATLAB, by itself, creates a new AUDIO_SIGN_PROCESSING_IN_MATLAB or raises the existing
%      singleton*.
%
%      H = AUDIO_SIGN_PROCESSING_IN_MATLAB returns the handle to a new AUDIO_SIGN_PROCESSING_IN_MATLAB or the handle to
%      the existing singleton*.
%
%      AUDIO_SIGN_PROCESSING_IN_MATLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUDIO_SIGN_PROCESSING_IN_MATLAB.M with the given input arguments.
%
%      AUDIO_SIGN_PROCESSING_IN_MATLAB('Property','Value',...) creates a new AUDIO_SIGN_PROCESSING_IN_MATLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Audio_Sign_Processing_in_matlab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Audio_Sign_Processing_in_matlab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Audio_Sign_Processing_in_matlab

% Last Modified by GUIDE v2.5 20-Oct-2018 20:47:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Audio_Sign_Processing_in_matlab_OpeningFcn, ...
                   'gui_OutputFcn',  @Audio_Sign_Processing_in_matlab_OutputFcn, ...
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


% --- Executes just before Audio_Sign_Processing_in_matlab is made visible.
function Audio_Sign_Processing_in_matlab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Audio_Sign_Processing_in_matlab (see VARARGIN)

% Choose default command line output for Audio_Sign_Processing_in_matlab
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Audio_Sign_Processing_in_matlab wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Audio_Sign_Processing_in_matlab_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Melody.
function Melody_Callback(hObject, eventdata, handles)
% hObject    handle to Melody (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clear,clc;
[data, fs] = audioread('orig_input.wav'); 
[rows colums] = size(data);
t = [0:1/fs:.20];
 
Al = sin(2*pi*220*t); 
A = sin(2*pi*440*t);
C = sin(2*pi*523.25*t);
D = sin(2*pi*587.33*t);
E = sin(2*pi*659.26*t);
H = sin(2*pi*261.6*t);
G = sin(2*pi*783.99*t);
Gl = sin(2*pi*196*t); 
B = sin(2*pi*493.88*t);
F = sin(2*pi*349.23*t);
Fl = sin(2*pi*174.61*t);
Rs = sin(2*pi*000*t); 
Dl = sin(2*pi*293.66*t); 
Gs = sin(2*pi*392*t); 
Bl = sin(2*pi*246.94*t);
El = sin(2*pi*164.81*t); 
As = sin(2*pi*466.16*t); 
 
AlA = (Al + A);
DAl = (D + Al);
CH = (C + H);
GGl=G+Gl;
FFl = (F + Fl);

y=[As,D,C,B,A,E,A,C,El,As,Bl,F,FFl,AlA,A,A,DAl,G,CH,GGl,G,B,C,D,Dl,A,As,C,F,A,F,B,A,FFl,C,A,D,DAl,As,A,B,G,A,Rs,C,G,B,CH,F,A,F,C,AlA,A,C,A,D,A,E,DAl,CH,C,E,C,G,C,E,CH,GGl,G,B,G,C,G,D,A,B,A,B]

data_new = y(1:length(data));

for i = 1:colums
    for j = 1:rows
        data_new(j+i) = data(j,i) + y(i+j);
        
    end
end

audiowrite('melody.wav', data_new, fs);
[z,Fs]=audioread('melody.wav');

soundsc(z,Fs);

 
% --- Executes on button press in FFT.
function FFT_Callback(hObject, eventdata, handles)
% hObject    handle to FFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
[y,fs] = audioread('melody.wav')
Y = fft(y);
plot(abs(Y))

N = 1024 % number of FFT points
transform = fft(y,N)/N;
magTransform = abs(transform);

faxis = linspace(-fs/2,fs/2,N);
plot(faxis,fftshift(magTransform));
xlabel('Frequency (Hz)')

% view frequency content up to half the sampling rate:
axis([0 length(faxis)/2, 0 max(magTransform)]) 

% change the tick labels of the graph from scientific notation to floating point: 
xt = get(gca,'XTick');  
set(gca,'XTickLabel', sprintf('%.0f|',xt))
% --- Executes on button press in spectogram.
function spectogram_Callback(hObject, eventdata, handles)
% hObject    handle to spectogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2);
[y,fs] = audioread('melody.wav')
win = 128 % window length in samples
% number of samples between overlapping windows:
hop = win/2            

nfft = win % width of each frequency bin 
spectrogram(y,win,hop,nfft,fs,'yaxis')

% change the tick labels of the graph from scientific notation to floating point: 
yt = get(gca,'YTick');  
set(gca,'YTickLabel', sprintf('%.0f|',yt))

% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close
