
function [] = analyze_simulation(pathResults)
  
  load(pathResults)
  
  offset_fr = 0;
  dfr = 0.2;
  i = 0;
  fr_lim = max(analysis.fr);
  analysis.res = zeros(ceil(fr_lim/dfr),5);
  
  while offset_fr < fr_lim
    i = i+1;
    min_fr = offset_fr;
    max_fr = offset_fr + dfr;
    
    idx = analysis.fr >= min_fr & analysis.fr < max_fr;
    analysis.res(i,:) = sum(analysis.accuracy(idx,:),1)/sum(idx);
    
    offset_fr = offset_fr + dfr;
  end
  
  figure
  subplot(2,1,1)
  histogram(analysis.fr,0:dfr:fr_lim)
  xlim([0,ceil(fr_lim)])
  
  subplot(2,1,2)
  hold on
  plot(dfr/2:dfr:fr_lim+dfr/2,analysis.res(:,2),'k')
  plot(dfr/2:dfr:fr_lim+dfr/2,analysis.res(:,3),'r')
  plot(dfr/2:dfr:fr_lim+dfr/2,analysis.res(:,4),'y')
  xlim([0,ceil(fr_lim)])
  ylim([0,0.5])
end