function [phaseStatsStruct, timingStruct] = ...
          calcMovementSequenceDataSummary(movementSequence, ...
              dataIndex, dataSeries, dataScale,numberOfPhases, ...
              numberOfInterpolationPoints)

assert(numberOfPhases==2);            
            
assert(isempty(dataSeries)==0);


numberOfRepetitions = 0;
for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexReference))==0)
    numberOfRepetitions=numberOfRepetitions+1;
  end  
end

rawDataStruct(numberOfRepetitions) = struct('x',[],'y',[]);

basicStatsStruct = struct('min',[],'p05',[],'p25',[],'median',[],'p75',[],...
                          'p95',[],'max',[],'n',0, 'data', rawDataStruct);      
                        
timingStruct(numberOfPhases) = struct('phase',basicStatsStruct);                        
                        
phaseStatsStruct(numberOfPhases) = struct('phase',basicStatsStruct,...
                                          'start',basicStatsStruct,...
                                          'end'  ,basicStatsStruct);
                        
                

                                        
normTime = [0:(1/(numberOfInterpolationPoints-1)):1]';                                        
                                        

%Initialize structs that require it
for p=1:1:numberOfPhases
  phaseStatsStruct(p).phase.n = 0;
  phaseStatsStruct(p).start.n = 0;
  phaseStatsStruct(p).end.n   = 0;
  timingStruct(p).phase.n     = 0;
end



for p=1:1:numberOfPhases
  for r=1:1:numberOfRepetitions
    phaseStatsStruct(p).phase.data(r).x = ...
      zeros(numberOfInterpolationPoints,1);
    phaseStatsStruct(p).phase.data(r).y = ...
      zeros(numberOfInterpolationPoints,1);
    
    phaseStatsStruct(p).start.data(r).x = zeros(1,1);
    phaseStatsStruct(p).start.data(r).y = zeros(1,1);
    
    phaseStatsStruct(p).end.data(r).x   = zeros(1,1);
    phaseStatsStruct(p).end.data(r).y   = zeros(1,1);
    
    timingStruct(p).phase.data(r).x           = zeros(1,1);
    timingStruct(p).phase.data(r).y           = zeros(1,1);
  end
end




idx = 1;
data(numberOfPhases)      = struct('val',[]);
dataStart(numberOfPhases) = struct('val',[]);
dataEnd(numberOfPhases)   = struct('val',[]);
dataTime(numberOfPhases)  = struct('val',[]);


for z=1:1:length(movementSequence)
  if( sum(isnan(movementSequence(z).indexReference))==0)
 
    idx0 = movementSequence(z).indexStart;
    idx1 = movementSequence(z).indexReference;
    idx2 = movementSequence(z).indexEnd;
    
    
    for p=1:1:numberOfPhases

      switch p
        case 1
          idxA = idx0;
          idxB = idx1;
        case 2
          idxA = idx1;
          idxB = idx2;          
        otherwise assert(0);          
      end
       
      %Make sure all of the values in this segment are actual numbers.
      assert( sum( isnan(dataIndex(idxA:idxB,1))) == 0);
      assert( sum( isinf(dataIndex(idxA:idxB,1))) == 0);
      
      %Make sure that this array is not just constants: something should be
      %varying
      assert( std(dataIndex(idxA:idxB,1)) > sqrt(eps));
      
      
      duration = (dataIndex(idxB,1)-dataIndex(idxA,1));      
      normDataTime = [dataIndex(idxA:1:idxB,1)-dataIndex(idxA,1)]./duration;
      
      
      data(p).val      = [data(p).val;      dataSeries(idxA:idxB,1).*dataScale];
      dataStart(p).val = [dataStart(p).val; dataSeries(idxA,1).*dataScale];
      dataEnd(p).val   = [dataEnd(p).val;   dataSeries(idxB,1).*dataScale];
      dataTime(p).val  = [dataTime(p).val;  duration];
      
      normData = interp1(normDataTime,...
                         dataSeries(idxA:idxB,1).*dataScale,...
                         normTime);
                       
      phaseStatsStruct(p).phase.data(idx).x = normTime.*duration+dataIndex(idxA,1);
      phaseStatsStruct(p).phase.data(idx).y = normData;
      phaseStatsStruct(p).phase.n = phaseStatsStruct(p).phase.n + idxB-idxA+1;
      
      phaseStatsStruct(p).start.data(idx).x = phaseStatsStruct(p).phase.data(idx).x(1,1);
      phaseStatsStruct(p).start.data(idx).y = phaseStatsStruct(p).phase.data(idx).y(1,1);
      phaseStatsStruct(p).start.n = 1+phaseStatsStruct(p).start.n;

      phaseStatsStruct(p).end.data(idx).x   = phaseStatsStruct(p).phase.data(idx).x(end,1);
      phaseStatsStruct(p).end.data(idx).y   = phaseStatsStruct(p).phase.data(idx).y(end,1);
      phaseStatsStruct(p).end.n = 1+phaseStatsStruct(p).end.n;

      timingStruct(p).phase.data(idx).x           = 1;
      timingStruct(p).phase.data(idx).y           = duration;
      timingStruct(p).phase.n = 1+timingStruct(p).phase.n;
    end
    idx=idx+1;
  end

end


  

for p=1:1:numberOfPhases  
  if(phaseStatsStruct(p).phase.n > 0)
    phaseStatsStruct(p).phase = ...
      calcBasicStatistics(data(p).val,phaseStatsStruct(p).phase);
    phaseStatsStruct(p).start = ...
      calcBasicStatistics(dataStart(p).val,phaseStatsStruct(p).start);
    phaseStatsStruct(p).end   = ...
      calcBasicStatistics(dataEnd(p).val,phaseStatsStruct(p).end);
    timingStruct(p).phase = calcBasicStatistics(dataTime(p).val,timingStruct(p).phase);

  end
end
  
