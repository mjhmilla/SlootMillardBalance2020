function [phaseResults,startResults,endResults] = calcGroupComparison( ...
                        groupMetricDataA, groupMetricDataB)

phaseResults = struct('p',NaN,'h',NaN,'alpha',NaN);
startResults = struct('p',NaN,'h',NaN,'alpha',NaN);
endResults   = struct('p',NaN,'h',NaN,'alpha',NaN);
                      
if( isfield(groupMetricDataA,'phase') && isfield(groupMetricDataB,'phase'))        
  if(isfield(groupMetricDataA.phase,'median') && isfield(groupMetricDataB.phase,'median'))
    dataA = [];    
    for i=1:1:length(groupMetricDataA.phase.data)
      %for k=1:1:size(groupMetricDataA.phase.data(i).y,2)
        dataA = [dataA; mean(groupMetricDataA.phase.data(i).y,2)];
      %end
    end
    dataB = [];    
    for i=1:1:length(groupMetricDataB.phase.data)
      %for k=1:1:size(groupMetricDataB.phase.data(i).y,2)
        dataB = [dataB; mean(groupMetricDataB.phase.data(i).y,2)];
      %end
    end
    [p,h] = ranksum(dataA,dataB);
    phaseResults.p = p;
    phaseResults.h = h;
    phaseResults.alpha = 0.05;    
  end
end

if( isfield(groupMetricDataA,'start') && isfield(groupMetricDataB,'start'))        
  if(isfield(groupMetricDataA.start,'median') && isfield(groupMetricDataB.start,'median'))
    dataA = [];    
    for i=1:1:length(groupMetricDataA.start.data)
      %for k=1:1:size(groupMetricDataA.start.data(i).y,2)
        dataA = [dataA; mean(groupMetricDataA.start.data(i).y,2)];
      %end
    end
    dataB = [];    
    for i=1:1:length(groupMetricDataB.start.data)
      %for k=1:1:size(groupMetricDataB.start.data(i).y,2)
        dataB = [dataB; mean(groupMetricDataB.start.data(i).y,2)];
      %end
    end
    [p,h] = ranksum(dataA,dataB);
    startResults.p = p;
    startResults.h = h;
    startResults.alpha = 0.05;    
  end
end

if( isfield(groupMetricDataA,'end') && isfield(groupMetricDataB,'end'))        
  if(isfield(groupMetricDataA.end,'median') && isfield(groupMetricDataB.end,'median'))
    dataA = [];    
    for i=1:1:length(groupMetricDataA.end.data)
      %for k=1:1:size(groupMetricDataA.end.data(i).y,2)
        dataA = [dataA; mean(groupMetricDataA.end.data(i).y,2)];
      %end
    end
    dataB = [];    
    for i=1:1:length(groupMetricDataB.end.data)
      %for k=1:1:size(groupMetricDataB.end.data(i).y,2)
        dataB = [dataB;mean(groupMetricDataB.end.data(i).y,2)];
      %end
    end
    [p,h] = ranksum(dataA,dataB);
    endResults.p = p;
    endResults.h = h;
    endResults.alpha = 0.05;    
  end
end