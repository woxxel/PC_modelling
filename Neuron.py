##% inputs

## p      - probability of being a place cell
## PC_dist  - distribution of PC locations (default: 'uniform')

## activity - structure with all kind of information about the firing activity
##      offset      - height and variance
##      sigma       - width and variance of the firing field
##      amplitude   - activity amplitude and variance of amplitude
##      num_fields
##      
## track    - structure with some information about the track
##      length      - length of the track
##      nbin        - number of bins

import numpy as np
import matplotlib.pyplot as plt
import math, random, time
import multiprocessing as mp

class Neuron:
  
  def __init__(self,nCells,track,activity,bh,PC_prob):
    
    self.PC_status = np.zeros(nCells).astype('bool')
    self.options = {'PC_dist':'uniform'}
    self.field = []
    self.bh = bh
    self.nCells = nCells
    
    self.activity = np.zeros((nCells,len(self.bh['time_raw'])))
    self.rate_pp = []
    self.x_max = 100
    
    for n in range(nCells):
      
      self.field.append({'nModes':np.NaN,
                         'theta':np.NaN,
                         'sigma':np.NaN,
                         'A':0,
                         'A_0':-1,
                         'activation':np.NaN})
      
      
      #self.PC_status[n] = (bool(activity['nModes']) & (random.random()<PC_prob))
      
      #if self.PC_status[n]:
      self.field[n] = self.set_PC_field(track,activity,PC_prob)
      if self.field[n]['nModes'] > 0:
        self.PC_status[n] = True
      #else:
        #self.field[n]['nModes'] = 0;
        #self.field[n]['theta'] = np.NaN;
        #self.field[n]['sigma'] = np.NaN;
        #self.field[n]['A'] = 0;
      
      #self.rate_pp.append(self.set_TC_fun(self.field[n]))
    
    
  def set_TC_fun(self,field):
    def TC_fun(x):
      
      x_max = 100
      TC = np.ones(len(x))*field['A_0']
      if field['nModes'] > 0:
        for j in [-1,0,1]:   ## loop, to have periodic boundary conditions
          TC += field['A']*np.exp(-(x-field['theta']+x_max*j)**2/(2*field['sigma']**2))
      return TC
    return TC_fun
      
  def set_PC_field(self,track,activity,PC_prob):
    
    if self.options['PC_dist'] == 'uniform':          ## uniform distribution
      nFields = random.randint(1,activity['nModes']) & (random.random()<PC_prob)
      field = {'nModes':nFields,
               'theta':np.random.rand(nFields)*track['nbin'] if nFields>0 else np.NaN,
               'sigma':np.zeros(nFields) if nFields>0 else np.NaN,
               'A':np.zeros(max(1,nFields)),
               'A_0':np.NaN,
               'activation':np.NaN}
      
      field['A_0'] = draft_para(activity['A_0'])
      if nFields >0:
        field['sigma'] = draft_para(activity['sigma'],nFields)
        field['A'] = draft_para(activity['A']*field['A_0'],nFields)
        field['activation'] = draft_para(activity['field_activation'],nFields)
      
      
    #elif self.options['PC_dist'] == 'salientgauss':     ## gaussian around salient locations provided
    
    #elif self.options['PC_dist'] == 'manual':           ## according to provided distribution
    return field
  
  def generate_activity_all(self,nP=0):
    
    t_start = time.time()
    if nP> 0:
      pool = mp.Pool(nP)
      self.activity = pool.map(self.generate_activity,range(self.nCells))
      #print(status)
      print('>>> all done. time passed: %6.4g secs <<<'%(time.time()-t_start))
    else:
      ### generate activity for neurons
      perc = 0.1;
      j = 0;
      for n in range(self.nCells):
        if n/self.nCells>=perc*j:
          print('>>> %3.0d%% done. time passed: %6.4g secs <<<'%(round(perc*j*100),time.time()-t_start))
          j=j+1;
        self.activity[n,:] = self.generate_activity(n)
    
  def generate_activity(self,n):
    # example of generating a 
    # nonhomogeneous poisson process on [0,T] with intensity function intens
    
    #print('processing neuron n=%d'%n)
    intens = self.set_TC_fun(self.field[n])#self.rate_pp[n];
    intens_pos = intens(self.bh['binpos_raw'])
    
    t = self.bh['time_raw']
    T_end = t[-1]
    dt = t[1]-t[0];
    max_intens = self.field[n]['A_0'] + self.field[n]['A'];   ## define rate as maximum rate of hompp
    
    u = np.random.rand(int(math.ceil(1.1*T_end*max_intens)));  ## generate a sufficient amount of random variables to cover the whole time
    
    t_AP = np.cumsum(-(1/max_intens)*np.log(u));       ## generate points of homogeneous pp (continuous time)
    t_AP = t_AP[t_AP<T_end];
    
    nAP=len(t_AP);
    idx_AP = np.zeros(nAP).astype('int')
    for (AP,i) in zip(t_AP,range(nAP)):
      idx_AP[i] = np.argmin(abs(AP-t))
      
    #idx_AP = np.floor(t_AP/dt).astype('int')
    #idx_AP = idx_AP[bh['longrunperiod'][idx_AP]]            ## kick out activity during non-running
    T_AP = t[idx_AP];                              ## points of homogeneous pp (discrete time)
    m = intens_pos[idx_AP]                         ## evaluates intensity function at homogeneous pp points
    
    #idx_AP_store =np.copy(idx_AP)
    trials_frame = np.hstack([0, np.where(np.diff(self.bh['binpos_raw'])<-10)[0]+1,len(self.bh['time_raw'])-1])
    
    if self.field[n]['nModes']>0:
      if self.field[n]['activation']<1:
        idx_keep = np.zeros(nAP).astype('bool')
        for i in range(self.bh['trials']['ct']):
          
          idx_trial_APs = np.where((T_AP>=t[trials_frame[i]]) & (T_AP<t[trials_frame[i+1]]))[0]
          if len(idx_trial_APs)==0:
            continue
          N_APs = len(idx_trial_APs)
          #print(idx_trial_APs)
          idx_first_trial_AP = idx_trial_APs[0]
          idx_last_trial_AP = idx_trial_APs[-1]
          
          #print([T_AP[idx_first_trial_AP],T_AP[idx_last_trial_AP]])
          
          if random.random() < self.field[n]['activation']:
            m_trial = m[idx_first_trial_AP:idx_last_trial_AP+1]
          else:
            m_trial = np.ones(N_APs)*m.min()
          
          idx_keep_trial = np.random.rand(N_APs)<(m_trial/max_intens)
          
          idx_keep[idx_first_trial_AP:idx_last_trial_AP+1] = idx_keep_trial#idx_trial_APs[idx_keep_trial])# = idx_AP[idx_keep]
        
      else:
        idx_keep = np.random.rand(nAP)<(m/max_intens)
      
      idx_AP = idx_AP[idx_keep];
      
      t_AP = t_AP[idx_keep];      ## filter out some points with prob according to ratio of nonhomo/homo pp
      T_AP = T_AP[idx_keep];
    nAP=len(t_AP);
    
    #print(np.ndtype(idx_AP))
    #print(idx_AP)
    for AP in idx_AP.astype('int'):
      self.activity[n,AP]+= random.gauss(0.2,0.05)
    #self.activity[n,idx_AP] = 1#np.random.gamma(2,0.5,nAP);     ## shape = k, scale = theta, mode = (k-1)*theta, mean = k*theta, var = k*theta^2
    
    plt_bool = False;
    if plt_bool:
      plt.figure()
      
      plt.subplot(222)
      x = np.linspace(0,100,101);
      plt.plot(x,intens(x),'k')
      
      plt.subplot(224)
      plt.hist(np.diff(T_AP),np.linspace(0,2,101))
      plt.plot(np.ones(2)/15,[0,10],'r')
      
      plt.subplot(321)
      plt.plot(self.bh['time_raw'],intens_pos,'k-')
      plt.scatter(T_AP,intens_pos[idx_AP],color='r',marker='o')
      
      plt.subplot(323)
      plt.plot(t,self.bh['binpos_raw'],'k')
      plt.scatter(T_AP,self.bh['binpos_raw'][idx_AP],color='r',marker='o')
      
      plt.subplot(325)
      plt.plot(t,self.activity[n,:],'r-')
      
      plt.show()
      plt.close('all')
    return self.activity[n,:]
    
def draft_para(in_range,n=1):
  return in_range[0] + (in_range[1] - in_range[0])*np.random.rand(n)
