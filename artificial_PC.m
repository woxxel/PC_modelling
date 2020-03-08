


function [neurons, analysis, PC_fields] = artificial_PC(nCells,s,pathResults)
  
  %  neurons = zeros(100,1);
  pathMouse='/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/M884';
  mouse='884';
  para = set_paras([],[],mouse);

  %  pathMouse = sprintf('%s%d',basePath,mouse);
  pathSession = pathcat(pathMouse,sprintf('Session%02d',s));
  pathBH = dir(pathcat(pathSession,'*aligned.mat'));
  pathBH = pathcat(pathSession,pathBH.name);
  %  if exist(pathSave,'file') && loadData
  %    disp(sprintf('Behaviour file already exists. Loading from %s...',pathSave))
  bh = load(pathBH);
  bh = bh.alignedData.resampled;
  disp('loading done')

  activity = struct('offset',[0.5,0.2],...
                    'width',[2,1],...
                    'amplitude',[10,3],...
                    'num_fields',1,...
                    'T',length(bh.time));

  track = struct('length',100,...
                'nbin',100);
                
                
  disp(sprintf('Initializing %d neurons...',nCells))
  tic
  %%% initialize neurons
  neurons = Neuron(nCells,track,activity,1.0,'uniform');
  
  
  %% construct poisson-spike trains
  dt_arr = [0 diff(bh.time)];
  T_end = bh.time(end);
  
  %%% based on rate of bin of last spike
  perc = 0.1;
  j = 0;
  
  for n = 1:nCells
%      neurons(n).firingmap
    if n./nCells>=perc*j
      disp(sprintf('>>> %3.0d%% done. time passed: %6.4g secs <<<',round(perc*j*100),toc))
      j=j+1;
    end
    
    neurons(n).activity = generate_activity(neurons(n),bh.binpos,bh.time);
%      if neurons(n).status.PC
%        field = neurons(n).firingmap >= mean(neurons(n).firingmap);
%        idx_field = find(field);
%        idx_start = idx_field(1);
%        idx_end = idx_field(end);
%        if idx_start > idx_end
%          tmp_idx_start = idx_start;
%          idx_start = idx_end;
%          idx_end = tmp_idx_start;
%          in_reward = bh.binpos > idx_start | bh.binpos < idx_end;
%        else
%          in_reward = bh.binpos > idx_start & bh.binpos < idx_end;
%        end
%        idx_in_reward = find(in_reward);
%        d_in_reward = diff(idx_in_reward);
%        idx_entry = find(d_in_reward > 1);
%        idx_binpos = idx_in_reward(idx_entry);
%      end
%      
%      T = 0;
%      i = 1;
%      while T < T_end
%        [~,idx] = min(abs(T-bh.time));
%        bin = bh.binpos(idx);
%        rate = neurons(n).firingmap(bin);
%        
%        dt = nextTime(rate);
%        [~,idx_next] = min(abs(T+dt - bh.time));
%        
%        if neurons(n).status.PC && (i < length(idx_binpos) && idx_next > idx_binpos(i))
%          T = bh.time(idx_binpos(i));
%          i = i + 1;
%        else
%          T = T + dt;
%          neurons(n).activity(idx) = neurons(n).activity(idx) + 1;
%        end
%        
%      end
  end
%    toc
  
  %%% based on dynamically adjusted firing rate (how??)
%    for i = 1:length(bh.time)
%      t = bh.time(i);
%      dt = dt_arr(i);
%      bin = bh.binpos(i);
%      p_arr = rand(nCells,1);
%      for n=1:nCells
%        neurons(n).activity(i) = floor(neurons(n).firingmap(bin)*dt*p_arr(n));
%      end
%    end
  
  
  disp(sprintf('Analyzing %d neurons...',nCells))
  tic
  
  plt = false;
  analysis = struct('accuracy',false(nCells,5),'fr',zeros(nCells,1));%,'MI',struct('frac',NaN,'value',NaN,'dist',zeros(1,para.repnum)*NaN,'binMI',zeros(1,para.nbin)*NaN));
  %%% 1: pp = positive-positive: correctly detected PC
  %%% 2: ppf = positive-positive-false: correctly detected PC, but in wrong position
  %%% 3: pn = positive-negative: falsely detected PC
  %%% 4: np = negative-positive: not-detected PC
  %%% 5: nn = negative-negative: correctly not-detected PC
  
  PC_fields(nCells) = build_PC_fields(para);
%    para
  
  i = 0;
%    accuracy = struct('pp',zeros(nCells,1),'ppf',zeros(nCells,1),'pn',zeros(nCells,1),'np',zeros(nCells,1),'nn',zeros(nCells,1));
  for n = 1:nCells
    
    if n./nCells>=perc*i
      disp(sprintf('>>> %3.0d%% done. time passed: %6.4g secs <<<',round(perc*i*100),toc))
      i=i+1;
    end
    
    PC_fields(n) = anaPC(neurons(n).activity,bh,para,'spikes',0.3);
    analysis.fr(n) = PC_fields(n).mean_fr;
    
    if PC_fields(n).status && neurons(n).status.PC
      if (PC_fields(n).fields.num == neurons(n).field.num) && (abs(mod(PC_fields(n).fields.center - neurons(n).field.center + para.nbin/2,para.nbin)-para.nbin/2) < 10)
        analysis.accuracy(n,1) = true;
      else
        analysis.accuracy(n,2) = true;
      end
    end
    if PC_fields(n).status && ~neurons(n).status.PC
      analysis.accuracy(n,3) = true;
%        figure
%        plot_ax_activity_artificial(gca,PC_fields(n),neurons(n),bh,[])
%        waitforbuttonpress
%        close all
    elseif ~PC_fields(n).status && neurons(n).status.PC
      analysis.accuracy(n,4) = true;
%        figure
%        plot_ax_activity_artificial(gca,PC_fields(n),neurons(n),bh,[])
%        waitforbuttonpress
%        close all
    elseif ~PC_fields(n).status && ~neurons(n).status.PC
      analysis.accuracy(n,5) = true;
%        plt = false;
    end
    
    if plt
      close all
      
      figure('position',[1000 500 1500 1200])
      fields = plot_ax_activity_artificial(gca,PC_fields(n),neurons(n),bh);
      
      waitforbuttonpress
    end
    
  end
  
  if nargin > 2
    %% saving results
    pathResults = pathcat(pathResults,sprintf('results_nC%d_s%d.mat',nCells,s))
    save(pathResults,'analysis','PC_fields','neurons','bh','-v7.3')
    
    analyze_simulation(pathResults)
  end
  
  
  
end


%%% now, run MI on this and see if PCs are detected

%% do analysis: which PCs are best detected (width, amplitude, firingrate, ...), number of fields, etc
%% enable multiple fields
%% baseline shouldn't be monotonous, rather random distributed?!


function dt = nextTime(rateParameter)
  dt = -log(1 - rand()) / rateParameter;
end


function [t_out] = generate_activity(neuron,pos,T)
  % example of generating a 
  % nonhomogeneousl poisson process on [0,T] with intensity function intens
  
  intens = neuron.rate_pp;
%    intens(linspace(0,60,7))
  
  t = linspace(T(1),T(end),length(T));
  dt = t(2)-t(1);
  
  max_intens = neuron.field.baseline + neuron.field.ampl;   %% define rate as maximum rate of hompp
  
  u = rand(1,ceil(1.1*T(end)*max_intens));  %% generate a sufficient amount of random variables to cover the whole time
  
  y = cumsum(-(1/max_intens)*log(u));       %% generate points of homogeneous pp
  y = y(y<T(end));                 %% select those points less than T
  n=length(y);
  
  y_idx = max(1,round(y/dt));
  t_y = t(y_idx);
  m = intens(pos(y_idx));              %% evaluates intensity function at homogeneous pp points
  
  kick = rand(1,n)<m/max_intens;
  y_idx = y_idx(kick);
  y = y(kick);      %% filter out some points with prob according to ratio of nonhomo/homo pp
  t_y = t_y(kick);
  
  t_out = zeros(1,length(T));
  t_out(y_idx) = 1;
  
  plt = false;
  if plt
    figure
    subplot(3,2,1)
    hold on
    plot(t,intens(pos),'k')
    scatter(t_y,intens(pos(y_idx)),'ro','filled')
    
    subplot(3,2,3)
    hold on
    plot(t,pos,'k')
    scatter(t_y,pos(y_idx),'ro','filled')
  %    hist(y,10)

    subplot(3,2,5)
    plot(t,t-T,'k')
    
    subplot(1,2,2)
    x = linspace(0,80,100);
    plot(x,intens(x),'k')
    waitforbuttonpress
    close all
  end
  
end