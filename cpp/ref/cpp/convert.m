%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Make Image from "*.txt" files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;

%read the .txt from the folder
load('output.txt');

%change .txt into image and save as .jpg
str='output';

% eval:use string as statement
imshow(eval(str)/10,'border','tight','initialmagnification','fit');

% 300,100 is the ordination, 500,500 is the height and width
set (gcf,'Position',[300,100,500,500]); 
axis normal; % normal is full size

MyColorMap=[
   105/255, 105/255, 105/255	%0_silver
   255/255, 255/255,   0/255	%1_yellow
   255/255,   0/255, 255/255	%2_orange
   255/255, 140/255,   0/255	%3_orange
    16/255, 180/255,  57/255	%4_orange
   255/255, 236/255,  62/255	%5_orange
    16/255, 180/255,  57/255	%6_orange
   255/255, 236/255,  62/255	%7_orange
    32/255, 178/255, 170/255	%8_light blue
   255/255,   0/255,   0/255  %9_dark blue
];
colormap(MyColorMap);

set(gcf,'PaperPositionMode','auto');
print('-djpeg','output');
