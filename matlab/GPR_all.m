%%
clc; clear; close all;
display('Extracting the data...');
   
minLat = 40.598850;
maxLat = 40.810476;
minLong = -112.001349;
maxLong =  -111.713403;

[tr_lat,tr_long,tr_time,tr_pm25,IDs,models] = readQueryFile('data/queriedData_july4th_events.csv');

tFreq = 1;
% Downsampling the data
% lat_tr  = lat_tr(:,1:6:end);
% long_tr = long_tr(:,1:6:end);
% tr_pm25  = tr_pm25(16:end,:);
% tr_time = tr_time(16:end);
% tr_pm25  = tr_pm25(1:24,:);
% tr_time = tr_time(1:24);
tr_pm25  = tr_pm25(130:140,:);
tr_time = tr_time(130:140);
% nt = size(PMS_tr,2);
nt  = length(tr_time);
nID = length(tr_lat);

%preprocess sensor data
display('Preprocessing the data...');

tr_elevs = longLat2Elevs(tr_long,tr_lat);

tr_elevs  = repmat(tr_elevs,1,nt);
tr_elevs = tr_elevs(:);
tr_time = repmat(tr_time',nID,1);
tr_time = tr_time(:);

%Building the grid of the query points
% nPtsMin = 20;
% stepSize = min((maxLong-minLong),(maxLat-minLat))/(nPtsMin-1);
minT = min(tr_time);
maxT = max(tr_time);
t=minT:0.5:maxT;
nLats = 10;
nLongs = 16;
lats  = linspace(minLat,maxLat,nLats);
longs = linspace(minLong,maxLong,nLongs);
% longs = minLong:stepSize:maxLong;
% lats = minLat:stepSize:maxLat;
% nLongs = length(longs);
% nLats = length(lats);
nts = length(t);

% ***NOTE THAT WE HAVE TO CONVERT ALL THE LONGS
%    AND LATS TOGETHER SO THEY SCALE PROPERLY
[Q_Long,Q_Lat,Q_T]=meshgrid(longs,lats,t);
[xh,xv] = longLat2Meter([reshape(Q_Long(:,:,1),nLongs*nLats,1);tr_long],[reshape(Q_Lat(:,:,1),nLongs*nLats,1);tr_lat]);
% [Q_xh,Q_xv] = longLat2Meter(reshape(Q_Long(:,:,1),nLongs*nLats,1),reshape(Q_Lat(:,:,1),nLongs*nLats,1));
xh = xh/1000;
xv = xv/1000;
Q_xh = xh(1:nLongs*nLats);
Q_xv = xv(1:nLongs*nLats);
tr_xh = xh(nLongs*nLats+1:end);
tr_xv = xv(nLongs*nLats+1:end);

Q_El = longLat2Elevs(reshape(Q_Long(:,:,1),nLongs*nLats,1),reshape(Q_Lat(:,:,1),nLongs*nLats,1));

Q_T =reshape(Q_T,nLongs*nLats*nts,1);
Q_xh = repmat(Q_xh,nts,1);
Q_xv = repmat(Q_xv,nts,1);
Q_El = repmat(Q_El,nts,1);
tr_xv  = repmat(tr_xv,1,nt);
tr_xv = tr_xv(:);
tr_xh = repmat(tr_xh,1,nt);
tr_xh = tr_xh(:);

tr_X=[tr_xh,tr_xv,tr_time];
% Q_X=[Q_xh,Q_xv,Q_El,Q_T];
Q_X=[Q_xh,Q_xv,Q_T;tr_X];
% Q_X=tr_X([27:nID:end,22:nID:end,45:nID:end,25:nID:end,62:nID:end],:);
%Cleaning the Nan values from the measurments
tr_pm25 = tr_pm25';
tr_X = tr_X(~isnan(tr_pm25),:);
tr_pm25Vec = tr_pm25(~isnan(tr_pm25));
clear lats longs Q_xh Q_xv tr_elevs tr_time tr_xv tr_xh xh xv
clear stepSize nt nPtsMin nID maxT minT minLat maxLat minLong maxLong
%% Applying the regression
% L0 = [1.5, 0.05, 4]; L=[0.1791,0.1791,3.0691,1.4822];sigmaF = 11.2996;
% inversion data
% L0 =[1,0.03,1];L = [0.4184, 0.4184, 0.2719, 1.0000];sigmaF = 9.0622
% inversion with no elev optimize
% L0 = [1.5, 0.03, 4];L =[0.0228,0.0228,0.0300,3.1771]; sigmaF=13.8603;
% One day 8:00-16:00
% L0=[1.5, 0.03, 1]; L=[0.0306,0.0306,0.0300,1.3992];sigmaF=15.4503
% One day 00:00-08:00
% L0=[1.5, 0.03, 1]; L=[0.0066,0.0066,0.0300,1.8990];sigmaF=15.9017
   
sigmaF0 = 10;%std(tr_pm25Vec);
L0 = [4.3, 4];
sigmaN = 4.2;

optL       = false;
optSigmaF  = false;
optSigmaN  = false;
center     = true;
effOpt     = true;
isARD      = true;
isSpatIsot = true;
learnRate  = 1e-4;
tol        = 1e-4;
maxIt      = 400;
basisFnDeg = 1;
[yPred,yVar] = gpRegression(tr_X,tr_pm25Vec,Q_X,tr_X,tr_pm25Vec,...
  sigmaF0,optSigmaF,L0,optL,sigmaN,optSigmaN,basisFnDeg,isARD,isSpatIsot,learnRate,tol,maxIt,effOpt,center);

% while (i*10000<length(Q_X)+10000)
%   display(num2str(i));
%   [yPredPart,yVarPart] = gpRegression(tr_X,tr_pm25Vec,...
%     Q_X((i-1)*10000+1:min(i*10000,length(Q_X)),:),tr_X,tr_pm25Vec,...
%     sigmaF0,optSigmaF,L0,optL,sigmaN,optSigmaN,basisFnDeg,isARD,...
%     isSpatIsot,learnRate,tol,maxIt,effOpt,center);
%   yPred = [yPred;yPredPart];
%   yVar = [yVar;yVarPart];
%   i = i+1;
% end
%%
% Plotting the results
tr_YPRED = yPred(length(Q_T)+1:end,:);
tr_YVAR = yVar(length(Q_T)+1:end,:);
YPRED = yPred(1:length(Q_T),:);
YVAR  = yVar(1:length(Q_T),:);
tr_YPRED = reshape(tr_YPRED,size(tr_pm25));
tr_YVAR = reshape(tr_YVAR,size(tr_pm25));
YPRED = reshape(YPRED,nLats,nLongs,nts);
YVAR = reshape(YVAR,nLats,nLongs,nts);

%%
scrsize = get(0,'Screensize');
figure('Position',[10,10,scrsize(3)-100,scrsize(4)-100]);
for i=1:length(tr_lat)
  p = subplot(7,10,i);
  plot(t(1:1:end),tr_YPRED(i,6:6:end-5),'b-',t(1:1:end),tr_pm25(i,6:6:end-5),'r.','LineWidth',2,'MarkerSize',25);
  set (gca,'FontSize',12);
end  

%% interactive scatter plot of the data
scrsize = get(0,'Screensize');

figure('Position',[scrsize(3)/6,10,scrsize(3)*2/3,scrsize(4)-100]);
p2=subplot(2,2,2);
Q_ElMat = reshape(Q_El(1:nLats*nLongs,1),nLats,nLongs);
pcolor(p2,Q_Long(:,:,1),Q_Lat(:,:,1),Q_ElMat);
colorbar(p2,'FontSize',16,'FontWeight','bold')
shading(p2,'interp');
set(p2,'FontSize',16,'FontWeight','bold');
colorbar(p2,'FontSize',16,'FontWeight','bold');
xlabel(p2,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p2,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(p2,'Elevation map','FontSize',16,'FontWeight','bold')

minYP = min(min(YPRED(:)),min(min(tr_pm25(:,6:6:end-5))) );
maxYP = max(max(YPRED(:)),max(max(tr_pm25(:,6:6:end-5))) );

p1=subplot(2,2,1);
scatter(p1,tr_long,tr_lat,30,tr_pm25(:,6),'Filled');
set(p1,'FontSize',16,'FontWeight','bold');
colorbar(p1,'FontSize',16,'FontWeight','bold')
grid(p1,'on');
axis(p1,[Q_Long(1,1) Q_Long(1,end) Q_Lat(1,1) Q_Lat(end,1)])
title(p1,'Data set that we do regression on','FontSize',16,'FontWeight','bold')
xlabel(p1,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p1,'y(lat) [km]','FontSize',16,'FontWeight','bold');
caxis(p1,[minYP maxYP]);
txt1=uicontrol('Style','text','FontSize',16,'FontWeight','bold',...
        'Position',[1100/2 (1100+70)/2 200 25],'HorizontalAlignment','center',...
        'String','data t=0');

p3=subplot(2,2,3);
pcolor(p3,Q_Long(:,:,1),Q_Lat(:,:,1),YPRED(:,:,1));

shading(p3,'interp');
hold(p3,'on');
scatter(p3,tr_long,tr_lat,30,'o','filled','MarkerFaceColor','r');
set(p3,'FontSize',16,'FontWeight','bold');
colorbar(p3,'FontSize',16,'FontWeight','bold');
xlabel(p3,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p3,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(p3,'predicted Y','FontSize',16,'FontWeight','bold')
caxis(p3,[minYP maxYP]);


p4=subplot(2,2,4);
pcolor(p4,Q_Long(:,:,1),Q_Lat(:,:,1),sqrt(YVAR(:,:,1)));
minVP = min(sqrt(YVAR(:)));
maxVP = max(sqrt(YVAR(:)));
shading(p4,'interp');
hold(p4,'on')
scatter(p4,tr_long,tr_lat,30,'o','filled','MarkerFaceColor','r');
set(p4,'FontSize',16,'FontWeight','bold');
colorbar(p4,'FontSize',16,'FontWeight','bold')
xlabel(p4,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p4,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(p4,'uncertainty','FontSize',16,'FontWeight','bold')
caxis(p4,[minVP maxVP]);


txt2=uicontrol('Style','text','FontSize',16,'FontWeight','bold',...
        'Position',[1100/2 (1100-60)/2 200 25],'HorizontalAlignment',...
        'center','String','reg. t=0');

sld = uicontrol('Style', 'slider',...
  'Min',0,'Max',nts-1,'Value',0,'SliderStep',[1/(nts-1) 2/(nts-1)],...
  'Position', [(1100)/2 (1100)/2 200 30],...
  'Callback', {@updateAllPlot,p1,p3,p4,tr_long,tr_lat,...
               Q_Long(:,:,1),Q_Lat(:,:,1),t,tr_pm25(:,6:6:end-5),YPRED,YVAR,...
               txt1,txt2,minYP,maxYP,minVP,maxVP});


%%
dupInds = [];
for i=1:length(tr_lat)
  for j=1:i-1
    if (tr_lat(j)==tr_lat(i) && tr_long(j)==tr_long(i))
      dupInds = [dupInds;[i,j]];
    end
  end
end
%%
scrsize = get(0,'Screensize');
figure('Position',[10,10,scrsize(3)-100,scrsize(4)-100]);
for i=1:length(dupInds)
  subplot(4,6,i);
  plot(t(1:1:end),tr_YPRED(dupInds(i,1),:),'b-',t(1:1:end),...
    tr_pm25(dupInds(i,1),:),'r.',t(1:1:end),tr_pm25(dupInds(i,2),:),'g.','LineWidth',2,'MarkerSize',15);
  title({['sensor IDs: ',IDs{dupInds(i,1)},' and ', IDs{dupInds(i,2)}],...
    ['Latitude: ', num2str(tr_lat(dupInds(i,1)))],...
    ['Longitude: ', num2str(tr_long(dupInds(i,1)))]},'FontWeight','normal');
  set (gca,'FontSize',10);
end

