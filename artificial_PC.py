import numpy as np
import time, sys, os, random, h5py
import Neuron

sys.path.append('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Programs/PC_analysis')
from get_fields import *#PC_detect, prepare_PC_detect, build_PC_results#,get_behavior
from find_PCs import *
from utils import pickleData

def generate_artificial_data(basePath,mouse,s,nCells,nP,plt_bool=False,sv_bool=False):
  
  pathMouse = pathcat([basePath,mouse])
  print('impact of randomized spike height?! shouldnt be all same?!')
  dPC = detect_PC(basePath,mouse,s,nP,plt_bool,sv_bool)
  
  pathSession = pathcat([pathMouse,'Session%02d'%s])
  pathResults = pathcat([pathSession,'artificialData.mat'])
  
  print('loading done')

  activity = {'A_0':[0,10],
              'sigma':[1,8],
              'A':[0.25,4],
              'nModes':1,
              'T':len(dPC.dataBH['time']),
              'field_activation':[0,1]}         ## probability for each trial, that the field is being activated

  track = {'length':100,
           'nbin':dPC.para['nbin']};
  
  
  print('Initializing %d neurons...'%nCells)
  t_start = time.time()

  ### initialize neurons
  neurons = Neuron.Neuron(nCells,track,activity,dPC.dataBH,0.8);
  print(dPC.dataBH['trials'].keys())
  
  ### generate activity for neurons
  perc = 0.1;
  j = 0;
  print('introduce chance of activation per trial!')
  for n in range(nCells):
    if n/nCells>=perc*j:
      print('>>> %3.0d%% done. time passed: %6.4g secs <<<'%(round(perc*j*100),time.time()-t_start))
      j=j+1;
    neurons.generate_activity(n,dPC.dataBH);
  
  neuronData = {}
  neuronData['fields'] = neurons.field
  neuronData['S'] = neurons.activity
  
  pickleData(neuronData,pathcat([pathSession,'artificialData.py']))
  #return neuronData, dPC


#def analyze_artificial_data(basePath,mouse,s):
  print('Analyzing %d neurons...'%nCells)
  t_start = time.time()
  PC_fields = dPC.run_detection(neurons.activity,return_results=True);
  
  pickleData(PC_fields,pathcat([pathSession,'artificialData_analyzed.py']))
  
  return PC_fields, neurons


def analyze_simulation(PC_fields,neurons):
  
  ## characterize results:
  ##% 1: pp = positive-positive: correctly detected PC
  ##% 2: ppf = positive-positive-false: correctly detected PC, but in wrong position
  ##% 3: pn = positive-negative: falsely detected PC
  ##% 4: np = negative-positive: not-detected PC
  ##% 5: nn = negative-negative: correctly not-detected PC
  
  nCells = PC_fields['fields']['parameter'].shape[0]
  print('adapt find_PCs, such that it can take these kind of values')
  
  analysis = {'accuracy':np.zeros((nCells,5)).astype('bool'),
              'fr':np.zeros(nCells)}
  
  idx_MI = fdr_control(PC_fields['status']['MI_p_value'],0.1);
  idx_Bayes = PC_fields['status']['Bayes_factor'][:,0]-PC_fields['status']['Bayes_factor'][:,1] > 0;
  
  PC_status = idx_MI & idx_Bayes
  print(PC_status)
  print(neurons.PC_status)
  i = 0;
#    accuracy = struct('pp',zeros(nCells,1),'ppf',zeros(nCells,1),'pn',zeros(nCells,1),'np',zeros(nCells,1),'nn',zeros(nCells,1));
  
  for n in range(nCells):
    
    #if n/nCells>=perc*i:
      #print('>>> %3.0d%% done. time passed: %6.4g secs <<<'%(round(perc*i*100),time.time()-t_start))
      #i=i+1;
    
    #print('fields:')
    #print(PC_fields['fields']['parameter'][n,0,3,0])
    #print(neurons.field[n]['theta'])
    
    if PC_status[n] & neurons.field[n]['nModes'] > 0:#neurons.PC_status[n]:
      if (PC_fields['fields']['nModes'][n] == neurons.field[n]['nModes']) & (abs(PC_fields['fields']['parameter'][n,0,3,0] - neurons.field[n]['theta']) <10):
      
      
      #((abs(PC_fields['fields']['parameter'][n,0,3,0] - neurons.field[n]['theta'] + dPC.para['nbin']/2)%dPC.para['nbin']-dPC.para['nbin']/2) < 10):
        analysis['accuracy'][n,0] = True;
      else:
        analysis['accuracy'][n,1] = True;
      
    if PC_status[n] & ~(neurons.field[n]['nModes'] > 0):#neurons.PC_status[n]:
      analysis['accuracy'][n,2] = True;
#        figure
#        plot_ax_activity_artificial(gca,PC_fields(n),neurons(n),bh,[])
#        waitforbuttonpress
#        close all
    elif ~PC_status[n] & (neurons.field[n]['nModes'] > 0):#neurons.PC_status[n]:
      analysis['accuracy'][n,3] = True;
#        figure
#        plot_ax_activity_artificial(gca,PC_fields(n),neurons(n),bh,[])
#        waitforbuttonpress
#        close all
    elif ~PC_status[n] & ~(neurons.field[n]['nModes'] > 0):#neurons.PC_status[n]:
      analysis['accuracy'][n,4] = True;
#        plt = false;
    
    #if plt_bool:
      #close all
      
      #plt.figure('position',[1000 500 1500 1200])
      #fields = plot_ax_activity_artificial(gca,PC_fields[n],neurons(n),bh);
      
      #waitforbuttonpress
  
  #if nargin > 2:
    ### saving results
    #pathResults = pathcat(pathResults,sprintf('results_nC%d_s%d.mat',nCells,s))
    #save(pathResults,'analysis','PC_fields','neurons','bh','-v7.3')
    
    #analyze_simulation(pathResults)
  
  
  return neurons, analysis, PC_fields


nP = 8
nCells = 100000
PCs, neurons = generate_artificial_data("/media/wollex/Analyze_AS3/Data","879",3,nCells,nP)
PCs, neurons = generate_artificial_data("/media/wollex/Analyze_AS3/Data","879",7,nCells,nP)
PCs, neurons = generate_artificial_data("/media/wollex/Analyze_AS3/Data","879",10,nCells,nP)

##% now, run MI on this and see if PCs are detected

## do analysis: which PCs are best detected (width, amplitude, firingrate, ...), number of fields, etc
## enable multiple fields
## baseline shouldn't be monotonous, rather random distributed?!

def pathcat(strings):
  return '/'.join(strings)


def fdr_control(x,alpha):
  
  x[x==0.001] = 10**(-10)
  x_mask = ~np.isnan(x)
  N = x_mask.sum()
  FDR_thr = range(1,N+1)/N*alpha
  
  idxes = np.where(x_mask)[0];
  idx_sorted = np.argsort(x[x_mask]);
  x_sorted = x[x_mask][idx_sorted]
  
  FDR_arr = x_sorted<FDR_thr;
  idx_cut = np.where(FDR_arr==False)[0][0]
  FDR_arr[idx_cut:] = False
  
  classified = np.zeros(N).astype('bool')
  classified[idxes[idx_sorted[FDR_arr]]] = True
  return classified
  
  return 