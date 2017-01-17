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
maxData = 10;   

proj = 'debug-20170117-134008';
files = dir(proj);

findIndex = @(y) str2double(y{1});
indices = arrayfun(@(x) findIndex(strsplit(x.name,'-')), files);
[~,I] = sort(indices);
files = files(I(1:end-2));

count = 0;
for file = files'
   
   cla;
   
   disp(file);
   
   %read the .txt from the folder
   data = load([proj,'/',file.name]);
   
   % eval:use string as statement
   imshow(data/maxData,'border','tight','initialmagnification','fit');
   
   % 300,100 is the ordination, 500,500 is the height and width
   set (gcf,'Position',[300,100,500,500]);
   axis normal; % normal is full size
   
   colormap(MyColorMap);
   
   pause(0.01);
   
   set(gcf,'PaperPositionMode','auto');
   print('-djpeg',[proj,'/',num2str(count),'.jpg']);
   count = count + 1;
   
end