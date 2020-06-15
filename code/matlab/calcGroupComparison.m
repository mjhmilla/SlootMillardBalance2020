function [phaseResults,startResults,endResults] = calcGroupComparison( ...
                        groupMetricDataA, groupMetricDataB)

phaseResults = struct('p',NaN,'h',NaN,'alpha',NaN);
startResults = struct('p',NaN,'h',NaN,'alpha',NaN);
endResults   = struct('p',NaN,'h',NaN,'alpha',NaN);
                      
if( isfield(groupMetricDataA,'phase') && isfield(groupMetricDataB,'phase'))        
  if(isfield(groupMetricDataA.phase,'median') && isfield(groupMetricDataB.phase,'median'))
    dataA = [];    
    for i=1:1:length(groupMetricDataA.phase.data)
      dataA = [dataA;groupMetricDataA.phase.data(i).y];
    end
    dataB = [];    
    for i=1:1:length(groupMetricDataB.phase.data)
      dataB = [dataB;groupMetricDataB.phase.data(i).y];
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
      dataA = [dataA;groupMetricDataA.start.data(i).y];
    end
    dataB = [];    
    for i=1:1:length(groupMetricDataB.start.data)
      dataB = [dataB;groupMetricDataB.start.data(i).y];
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
      dataA = [dataA;groupMetricDataA.end.data(i).y];
    end
    dataB = [];    
    for i=1:1:length(groupMetricDataB.end.data)
      dataB = [dataB;groupMetricDataB.end.data(i).y];
    end
    [p,h] = ranksum(dataA,dataB);
    endResults.p = p;
    endResults.h = h;
    endResults.alpha = 0.05;    
  end
end