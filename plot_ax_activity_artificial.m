
function plot_ax_activity_artificial(ax,PC_fields,neurons,bh,plot_options)
  
  para = set_paras([],[]);
                
  hold(ax,'on')
  
%    fields = struct('occ',zeros(para.nbin,1));
%    extent = struct('start_p',0,'end_p',0);
  
%    if PC_fields.status > 0
    for f = 1:PC_fields.fields.num
        
        if PC_fields.fields.frac_active(f) > 0.3
            col = [1, 0.8, 0.8];
        else
            col = [1, 0.95, 0.95];
        end
        idx = 1:para.nbin;
        
        if PC_fields.fields.extents(f,2) > PC_fields.fields.extents(f,1)
          idx = idx > PC_fields.fields.extents(f,1) & idx < PC_fields.fields.extents(f,2);
%            barh(ax,(PC_fields.fields.map(PC_fields.fields.extents(f,1):PC_fields.fields.extents(f,2))>0)*bh.time(end),1,'FaceColor',col,'EdgeColor','None')
        else
          idx = idx > PC_fields.fields.extents(f,1) | idx < PC_fields.fields.extents(f,2);
%            barh(ax,(PC_fields.fields.map(PC_fields.fields.extents(f,2):end)>0)*bh.time(end),1,'FaceColor',col,'EdgeColor','None')
%            barh(ax,(PC_fields.fields.map(1:PC_fields.fields.extents(f,1))>0)*bh.time(end),1,'FaceColor',col,'EdgeColor','None')
        end
        barh(ax,idx*bh.time(end),1,'FaceColor',col,'EdgeColor','None')
%          fields(f).occ = PC_fields.fields.map;
    end
%    end
  
  if neurons.status.PC
    barh(ax,bh.time(end)*(neurons.firingmap > mean(neurons.firingmap)),1,'FaceColor',[0.9,0.9,0.9],'EdgeColor','none');
  end
  
  barh(ax,-neurons.firingmap*(bh.time(end)/10)/max(PC_fields.firingmap),'k','FaceAlpha',1)
  barh(ax,-PC_fields.firingmap*(bh.time(end)/10)/max(PC_fields.firingmap),'r','FaceAlpha',0.9)
  
  plot(ax,bh.time,bh.position/para.binwidth)
%    neurons.activity = neurons.activity / (prctile(neurons.activity(neurons.activity>0),20)*2);
  
  if length(neurons.activity) > 0
    idx = neurons.activity & bh.longrunperiod;
    activity_lr = zeros(1,8989);
    activity_lr(idx) = neurons.activity(idx);
    scatter(ax,bh.time(find(activity_lr>0)),bh.position(activity_lr>0)/para.binwidth,3*activity_lr(activity_lr>0)+5,'r','fill')
    
    idx = neurons.activity & ~bh.longrunperiod;
    activity_nlr = zeros(1,8989);
    activity_nlr(idx) = neurons.activity(idx);
    scatter(ax,bh.time(find(activity_nlr>0)),bh.position(activity_nlr>0)/para.binwidth,3*activity_nlr(activity_nlr>0)+5,'k','fill')
  end
%    xlims = xlim;
%    xlim(ax,[xlims(1) 8989./15])
  hold(ax,'off')
  set(ax,'YAxisLocation','right')
  ylim(ax,[0 para.nbin])
  
  set(ax,'YTick',linspace(0,80,5))
  set(ax,'YTickLabels',linspace(0,100,5))
  xlabel(ax,'t [s]')
  ylabel(ax,'Location [cm]')
  
  if nargin == 5
    for i = 1:length(plot_options)
      switch plot_options{i}
        case 'no_xlabel'
          xlabel('')
        case 'no_xticks'
          xticks([])
        case 'no_ylabel'
          ylabel([])
        case 'no_yticks'
          yticks([])
      end
    end
  end
  
  
%    smoothed_MI = imgaussfilt(PC_fields.MI.binMI,para.sigma,'Padding','circular','FilterDomain','spatial');
%    periodic_MI = [smoothed_MI(end-(para.offset-1):end) smoothed_MI smoothed_MI(1:para.offset)];
%    MI_thr = 2*nanmean(smoothed_MI);
  
  smoothed_fr = imgaussfilt(PC_fields.firingmap,para.sigma,'Padding','circular','FilterDomain','spatial');
  periodic_fr = [smoothed_fr(end-(para.offset-1):end) smoothed_fr smoothed_fr(1:para.offset)];
  
  fr_baseline = prctile(smoothed_fr,50);
  fr_baseline = mean(smoothed_fr(smoothed_fr<=fr_baseline));
  
  fr_max = max(smoothed_fr);
  fr_thr = fr_baseline + 0.25*(fr_max-fr_baseline);
  
  fr_mean = mean(PC_fields.firingmap);
  fr_thr2 = 1.25*fr_mean;
  
%    disp('vars:')
%    sqrt(var(PC_fields.firingmap))
%    sqrt(var(PC_fields.MI.binMI))
  
%    fields = bwareaopen(periodic_fr > max(fr_thr,fr_thr2),4);
%    smoothed_field = fields(para.offset+1:end-para.offset).*smoothed_MI;
%    if plt
  figure('position',[100 100 600 400])
%    subplot(2,1,1)
  hold on
  bar(PC_fields.firingmap,'r')
  bar(smoothed_fr,'g','FaceAlpha',0.5)
  
  plot([0,80],[fr_baseline,fr_baseline],'k:')
  plot([0,80],[fr_thr,fr_thr],'k-')
  plot([0,80],[fr_max,fr_max],'k--')
  
  plot([0,80],[fr_mean,fr_mean],'k:','Linewidth',2)
  plot([0,80],[fr_thr2,fr_thr2],'k-','LineWidth',2)
  
  xlim([-10,90])
  
  
  fields = bwareaopen(periodic_fr > max(fr_thr,fr_thr2),4);
  
  
  [bw_comp,~] = bwlabel(fields);
      
  %% label twice appearing fields by the same label
  entries = mod(find(bw_comp)-1,para.nbin)+1;
  multi_idx = sum(entries==entries');
  for i = 1:length(entries)
      if multi_idx(i) > 1 && entries(i)<2*para.offset
          val = bw_comp(entries(i));
          old_val = bw_comp(entries(i)+para.nbin);
          if ~(val == old_val)
            bw_comp(find(bw_comp==old_val)) = val;
          end
      end
  end
  
%    bw_comp
%    unique(bw_comp)
  
  trial_thr = 0;
  mode_spikes = true;
  activity_lr = neurons.activity(bh.longrunperiod);                                               %% only activity from actual times
  spike_times = find(activity_lr);
  spikes = activity_lr(spike_times);
  mean_session = sum(spikes)/(sum(bh.longrunperiod)/para.f);
  j = 0;
  for i = unique(bw_comp)
    if i == 0
      continue
    end
    
    loc = find(bw_comp==i)-para.offset;
    loc_start = loc(1);
    loc_end = loc(end);
    
    if loc_start < 1 || loc_start > para.nbin
      loc_start = mod(loc_start-1,para.nbin)+1;
    end
    if loc_end < 1 || loc_end > para.nbin
      loc_end = mod(loc_end-1,para.nbin)+1;
    end
    
    %% check, whether cell is active in at least 40% of the trials
    active_in_field_ct = 0;
    for t = 1:PC_fields.trials.ct
      t_start = PC_fields.trials.frame(t)+1;
      t_end = PC_fields.trials.frame(t+1);
      pos = bh.binpos(t_start:t_end);
      longrun = bh.longrunperiod(t_start:t_end);
      act = neurons.activity(t_start:t_end);    % this should be in the field!
      
      if loc_start > loc_end
        mouse_in_field = (pos > loc_start | pos < loc_end);
      else
        mouse_in_field = (pos > loc_start & pos < loc_end);
      end
      
      act_field = act(mouse_in_field & longrun);
      if mode_spikes
        mean_trial = sum(act_field)/(sum(longrun(mouse_in_field))/para.f);
      else
        mean_trial = mean(act_field);
      end
      
      PC_fields.trials.PC_rate(t) = mean_trial;
      
      if mean_trial > 1.5*mean_session
        active_in_field_ct = active_in_field_ct + 1;
      end
    end
    
    frac_active = active_in_field_ct/PC_fields.trials.ct;
%          disp(sprintf('trials > average: %5.3g',frac_active))
    
    if active_in_field_ct > 1 && frac_active >= trial_thr
      j = j+1;
      
%            loc = loc + para.offset;
      idx = 1:para.nbin;
      if loc_start > loc_end
        idx_field = idx > loc_start | idx < loc_end;
        
        loc_start_shift = mod(loc_start+para.nbin/2-1,para.nbin)+1;
        loc_end_shift = mod(loc_end+para.nbin/2-1,para.nbin)+1;
        idx_field_shift = idx > loc_start_shift & idx < loc_end_shift;
        
        PC_fields.fields.width(j) = loc_end_shift - loc_start_shift + 1;
        PC_fields.fields.center(j) = mod(round(dot(idx(idx_field_shift),smoothed_fr(idx_field))/sum(smoothed_fr(idx_field)))-para.nbin/2-1,para.nbin)+1;
      else
        idx_field = idx > loc_start & idx < loc_end;
        
        PC_fields.fields.width(j) = loc_end - loc_start + 1;
        PC_fields.fields.center(j) = mod(round(dot(idx(idx_field),smoothed_fr(idx_field))/sum(smoothed_fr(idx_field)))-1,para.nbin)+1;
      end
      
      PC_fields.fields.frac_active(j) = active_in_field_ct/PC_fields.trials.ct;
      PC_fields.fields.extents(j,:) = [loc_start,loc_end];
%            PC_fields.fields.map(j,:) = zeros(1,para.nbin);
%            PC_fields.fields.map(j,idx_field) = smoothed_field(idx_field);
      
      PC_fields.status = true;
    end
  end
  
  PC_fields.fields.num = j;
%    PC_fields.fields.map = smoothed_field;
  
%    PC_fields
%    PC_fields.fields
%    PC_fields.MI
%    PC_fields.trials.PC_rate
  
  
      
  
  
%      waitforbuttonpress
%      close all
%    end
  
    
%    figure
%    histogram(PC_fields.MI.dist)
end