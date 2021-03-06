function capData = processCapturePoint(...
  mass,...
  wholeBodyData,...
  columnCenterOfMassPosition,...
  columnCenterOfMassVelocity,...
  gravityVector,...
  contactPlanes,...
  capFileName,...
  flag_loadFromFile)

if(flag_loadFromFile==0)
  npts = size(wholeBodyData,1);
  nPlanes = size(contactPlanes,1);  
  capData(nPlanes) = struct( 'r0F0', zeros(npts,3),...
                  'r0G0', zeros(npts,3),...  
                    'lu', zeros(npts,1),...
                    'vu', zeros(npts,1),...
                    'vk', zeros(npts,1),...
                     'u', zeros(npts,3),...
                     'n', zeros(npts,3),... 
                     'k', zeros(npts,3),...
                     'h', zeros(npts,1),...
                  'eorb', zeros(npts,1));
  for i=1:1:npts
    r0C0 = wholeBodyData(i,columnCenterOfMassPosition)';
    v0C0 = wholeBodyData(i,columnCenterOfMassVelocity)';
    
    
    for j=1:1:size(contactPlanes,1)
      capInfo = calcCapturePointInfo( mass,...
                                      r0C0,...
                                      v0C0,...
                                      contactPlanes(j,:)',...
                                      gravityVector); 

      capData(j).r0F0(i,:) = capInfo.r0F0';    
      capData(j).r0G0(i,:) = capInfo.r0G0';    
      capData(j).lu(i,1)   = capInfo.lu   ;  
      capData(j).vu(i,1)   = capInfo.vu   ;  
      capData(j).vk(i,1)   = capInfo.vk   ;  
      capData(j).u(i,:)    = capInfo.u'   ; 
      capData(j).n(i,:)    = capInfo.n'   ;    
      capData(j).k(i,:)    = capInfo.k'   ; 
      capData(j).h(i,1)    = capInfo.h    ; 
      capData(j).eorb(i,1) = capInfo.eorb ;                                        
    end
  end
  

  save(capFileName,...
      'capData');              

else
  load(capFileName);  
end