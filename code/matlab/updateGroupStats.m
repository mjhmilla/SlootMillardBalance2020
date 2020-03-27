function groupData = updateGroupStats(metricName,...
                                      groupMembers,...
                                      subjectData,...
                                      groupData,...
                                      flag_verbose)




n=0;




subFields = {'phase','start','end'};

for indexSubField = 1:1:length(subFields)
  dataStruct(length(groupMembers)) = struct('x',[],'y',[]);
  data = [];  
  for i = 1:1:length(groupMembers)
    indexSubject = groupMembers(i);  
    fieldName = subFields{indexSubField};
    
    if( isfield( subjectData(indexSubject).(metricName),fieldName)==1)            
      if(length( subjectData(indexSubject).(metricName).(fieldName).data) >0)
        

        %Go get the average data value for this subject
        avgX = zeros(size(subjectData(indexSubject).(metricName).(fieldName).data(1).x));
        avgY = zeros(size(subjectData(indexSubject).(metricName).(fieldName).data(1).y));
        for(k=1:1:length(subjectData(indexSubject).(metricName).(fieldName).data))
          avgX = avgX + subjectData(indexSubject).(metricName).(fieldName).data(k).x;      
          avgY = avgY + subjectData(indexSubject).(metricName).(fieldName).data(k).y;                     
        end
        avgX = avgX ./ length(subjectData(indexSubject).(metricName).(fieldName).data);
        avgY = avgY ./ length(subjectData(indexSubject).(metricName).(fieldName).data);

        %Store the average value as an entry into the data category of the group
        dataStruct(i).x = avgX;
        dataStruct(i).y = avgY;
        data = [data;avgY];
      end      
    end
  end
  
  groupData.(metricName).(fieldName).data = dataStruct;
  groupData.(metricName).(fieldName) = ...
    calcBasicStatistics( data, ...
      groupData.(metricName).(fieldName));

end

                  
if(flag_verbose==1)
 here=1;
end