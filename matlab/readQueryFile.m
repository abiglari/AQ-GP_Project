function [lats,longs,times,pm2_5,IDs,models] = readQueryFile(fileName, neglectMissing)

if (nargin<2)
  neglectMissing = false;
end
fid = fopen(fileName,'r');
rows = textscan(fid,'%s');
fclose(fid);
rows=rows{1};
nr = length(rows);

IDs = textscan(rows{1},'%q','delimiter', ',');
IDs = IDs{1}(3:end);
nIDs = length(IDs);

models = textscan(rows{2},'%q','delimiter', ',');
models = models{1}(3:end);

latsTxt = textscan(rows{3},'%q','delimiter', ',');
latsTxt = latsTxt{1}(3:end);

longsTxt = textscan(rows{4},'%q','delimiter', ',');
longsTxt = longsTxt{1}(3:end);

lats  = zeros(length(latsTxt),1);
longs = zeros(length(longsTxt),1);
remCols = [];
for i=1:length(latsTxt)
  lats(i)  = str2double(latsTxt{i});
  longs(i) = str2double(longsTxt{i});
  for j=1:i-1
    if lats(i)==lats(j) && longs(i)==longs(j)
      remCols = [remCols,i];
    end
  end
end
  
nt = nr-4;
pm2_5 = zeros(nt,nIDs);
times = zeros(nt,1);
for r=5:nr
  row = textscan(rows{r},'%q','delimiter', ',');
  date =  datetime(row{1}{1},'InputFormat','uuuu-MM-dd''T''HH:mm:ss''Z');
  if r==5
    startDate = date;
  end
  times(r-4) = (date - startDate) / hours(1);
  txts = row{1}(3:end);
  for c=1:length(txts)
    if (isempty(txts{c}) || str2double(txts{c})>300)
      pm2_5(r-4,c) = -1;
    else
      pm2_5(r-4,c) = str2double(txts{c});
    end
  end
end

% pm2_5(:,remCols) = [];
% models(remCols)  = [];
% IDs(remCols)     = [];
% lats(remCols)    = [];
% longs(remCols)   = [];
% nIDs = length(IDs);

j=1;
while j<=nIDs
  if sum(pm2_5(pm2_5(:,j)>0,j)==pm2_5(find(pm2_5(:,j)>0,1),j)) == sum(pm2_5(:,j)>0)
    pm2_5 = [pm2_5(:,1:j-1),pm2_5(:,j+1:end)];
    models = [models(1:j-1);models(j+1:end)];
    IDs = [IDs(1:j-1);IDs(j+1:end)];
    lats = [lats(1:j-1);lats(j+1:end)];
    longs = [longs(1:j-1);longs(j+1:end)];        
    nIDs = nIDs-1;
  else
    j = j+1;
  end
end

if (~neglectMissing)
  for j=1:nIDs
    if sum(pm2_5(:,j)<=0)~=0
      i =0;
      while i<nt
        i = i+1;
        if pm2_5(i,j)<=0 && i==1
          ind = find(pm2_5(:,j)>0,1);
          pm2_5(1:ind-1,j) = NaN;
          i = ind-1;
        elseif pm2_5(i,j)<=0
          z=i+1;
          while z<=nt && pm2_5(z,j)<=0
            z=z+1;
          end
          if (z<=nt && (z-i)<=6)
            while (i<z)
              pm2_5(i,j) = pm2_5(i-1,j) + (pm2_5(z,j)-pm2_5(i-1,j))/(z-(i-1));
              i=i+1;
            end
          else
            pm2_5(i:z-1,j) =NaN;
          end
          i=z;
        end
      end
    end
  end
end

for c=1:nIDs
  pm2_5(:,c) = calibrate(pm2_5(:,c),models{c});
end
end
