%%% inputs

%% p      - probability of being a place cell
%% PC_dist  - distribution of PC locations (default: 'uniform')

%% activity - structure with all kind of information about the firing activity
%%      offset      - height and variance
%%      width       - width and variance of width of the firing field
%%      amplitude   - activity amplitude and variance of amplitude
%%      num_fields
%%      
%% track    - structure with some information about the track
%%      length      - length of the track
%%      nbin        - number of bins

classdef Neuron < handle
  
  properties
    status;       %%% contains, whether the neuron is a place cell or not
    options;       %%% structure containing a few options
          %% PC_dist - density distribution for firing fields
          %% PC_shape - shape of firing field
    field;
%      firingmap;
    activity;
    rate_pp;
  end
  
  
  methods
    
    function obj = Neuron(nCells,track,activity,varargin)
      
      if nargin ~= 0
        
        for n = 1:nCells
          if nargin < 5
            obj(n).options.PC_dist = 'uniform'
          else
            obj(n).options.PC_dist = varargin{2};
          end
          
%            obj(n).firingmap = max(0,activity.offset(1)+(activity.offset(2)-activity.offset(1))*rand())*ones(1,track.nbin);
          
          obj(n).field.baseline = -1;
          while obj(n).field.baseline < 0
            obj(n).field.baseline = activity.offset(1)+(activity.offset(2)-activity.offset(1))*rand();
          end
          
%            obj(n).firingmap(1)
%            normrnd(activity.offset(1),activity.offset(2)))*ones(1,track.nbin);
          
          PC_status_tmp = activity.num_fields && rand()<varargin{1};
          obj(n).status = struct('PC',PC_status_tmp);
          
          if obj(n).status.PC
            obj(n).set_PC_field(track,activity)
          else
            obj(n).field.num = 0;
            obj(n).field.center = NaN;
            obj(n).field.width = NaN;
            obj(n).field.ampl = 0;
          end
          obj(n).activity = zeros(1,activity.T);
          
%            obj(n).field
          if ~obj(n).status.PC
            obj(n).rate_pp = @(x)(obj(n).field.baseline*ones(1,length(x)));
          elseif obj(n).field.num == 1
            obj(n).rate_pp = @(x)(obj(n).field.baseline*ones(1,length(x)) + ...
                                  obj(n).field.ampl(1)*exp(-(x-obj(n).field.center(1)).^2/(2*obj(n).field.width(1)^2)));
          elseif obj(n).field.num == 2
            %% dont have that anyway...
            obj(n).rate_pp = @(x)(obj(n).field.baseline*ones(1,length(x)) + ...
                                  obj(n).field.ampl(1)*exp(-(x-obj(n).field.center(1)).^2/(2*obj(n).field.width(1)^2)) + ...
                                  obj(n).field.ampl(2)*exp(-(x-obj(n).field.center(2)).^2/(2*obj(n).field.width(2)^2)));
          end
        end
      end
    end 
    
    
    function set_PC_field(obj,track,activity)
      
      switch obj.options.PC_dist
        case 'uniform'          %% uniform distribution
          
          nFields = randi(activity.num_fields);
          obj.field.num = nFields;
          obj.field.center = rand(nFields,1)*track.nbin;
          
          obj.field.width = zeros(nFields,1);
          while any(obj.field.width < 1)
            obj.field.width = normrnd(activity.width(1),activity.width(2),nFields,1);
          end
          
          obj.field.ampl = zeros(nFields,1);
          while any(obj.field.ampl < 1)
            obj.field.ampl = normrnd(activity.amplitude(1),activity.amplitude(2),nFields,1);
          end
%            obj.field
          
%            for f = 1:nFields
%              c = obj.field.center(f);
            
            
            %%% adjust for circular
%              field = normpdf(1:track.nbin,c-track.nbin,obj.field.width(i));
%              field = field + normpdf(1:track.nbin,c,obj.field.width(i));
%              field = field + normpdf(1:track.nbin,c+track.nbin,obj.field.width(i));
            
%              obj.firingmap = obj.firingmap + obj.field.ampl*field/max(field);
            
%            end
          
        case 'salientgauss'     %% gaussian around salient locations provided
        
        case 'manual'           %% according to provided distribution
      
      end
      
    end
    
  end
  
  
  
  
  
end