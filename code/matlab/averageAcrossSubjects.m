function summaryData = averageAcrossSubjects(group,...
                                      subjectSummaryData,...
                                      summaryDataStructName,...
                                      indexTrialType,...
                                      flag_verbose)

meanTrajectory
%n = length(group.index);    
n=0;
for indexSubject = 1:1:length(subjectSummaryData)
  data=subjectSummaryData(indexSubject).(summaryDataStructName)(:,indexTrialType);  
  if(sum(isnan(data))==0)
    if( any(group.index == indexSubject) )
      summaryData = summaryData  + data;
      n=n+1;
    end
  end
end

summaryData = summaryData./n;

if(flag_verbose==1)
  fprintf('  %s\n',summaryDataStructName);
  fprintf('  : %f min\n  : %f mean\n  : %f max\n  : %f 25\n  : %f 75\n  : %f SO\n  : %f ST\n',...
          summaryData(1,1),summaryData(3,1),summaryData(5,1),...
          summaryData(2,1),summaryData(4,1),...
          summaryData(6,1),summaryData(7,1));
end