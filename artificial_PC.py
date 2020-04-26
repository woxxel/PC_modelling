import numpy as np
import time, sys, os, random, h5py
import Neuron
import pickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

sys.path.append('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Programs/PC_analysis')
from get_fields import *#PC_detect, prepare_PC_detect, build_PC_results#,get_behavior
from find_PCs import *
from utils import pickleData, get_nPaths, extend_dict, pathcat, fdr_control

class artificialData:
  
  
  def __init__(self,basePath,mouse,s,nP,redo=False,plt_bool=False,sv_bool=False):
    
    self.mouse = mouse
    self.s = s
    self.basePath = basePath
    self.pathMouse = pathcat([basePath,mouse])
    self.pathSession = pathcat([self.pathMouse,'Session%02d'%s])
    self.pathNeurons = pathcat([self.pathSession,'artificialData.py'])
      
    self.dPC = detect_PC(basePath,mouse,s,nP,plt_bool,sv_bool)   ## obtain behavior data
    self.nP = nP
    self.neuronData = {}
    self.redo = redo
    
    
  def load_activity(self):
    
    f = open(self.pathNeurons,'rb')
    self.neuronData = pickle.load(f)
    f.close()
    
    nDat,pathData = get_nPaths(self.pathSession,'artificialData_analyzed_n')
    self.PC_fields = {}
    for i in range(nDat):
      print(pathData[i])
      f = open(pathData[i],'rb')
      ld_tmp = pickle.load(f)
      f.close()
      self.PC_fields = extend_dict(self.PC_fields,ld_tmp['fields']['parameter'].shape[0],ld_tmp)
      
  def generate_activity(self,nCells):
    
    print('impact of randomized spike height?! shouldnt be all same?!')
    
    activity = {'A_0':[0,10],
                'sigma':[0.5,10],
                'A':[0.1,5],
                'nModes':1,
                'T':len(self.dPC.dataBH['time']),
                'field_activation':[0.3,1]}         ## probability for each trial, that the field is being activated
    
    track = {'length':100,
            'nbin':self.dPC.para['nbin']};
    
    
    if self.redo | (not os.path.exists(self.pathNeurons)):
      print('Initializing %d neurons...'%nCells)
      
      ### initialize neurons
      neurons = Neuron.Neuron(nCells,track,activity,self.dPC.dataBH,0.8);
      neurons.generate_activity_all(self.nP);
      
      self.neuronData['fields'] = neurons.field
      self.neuronData['S'] = sp.sparse.csr_matrix(neurons.activity)
      
      pickleData(self.neuronData,self.pathNeurons,'save')
    else:
      self.neuronData = pickleData([],self.pathNeurons,'load')

  def analyze(self,batchSize=1000,rerun=False):
    
    nCells = self.neuronData['S'].shape[0]
    print('Analyzing %d neurons...'%nCells)
    t_start = time.time()
    
    nBatches = int(np.ceil(nCells/batchSize))
    self.PC_fields = {}
    for i in range(nBatches):
      print('\t\t\t ------ mouse %s --- session %d ------ %d / %d neurons processed\t ------ \t time passed: %7.2fs'%(self.mouse,self.s,i*batchSize,nCells,time.time()-t_start))
      pathSave = pathcat([self.pathSession,'artificialData_analyzed_n%02d.py'%(i+1)])
      if self.redo | rerun | (not os.path.exists(pathSave)):
        
        self.dPC.para['svname_art'] = pathSave
        PC_fields_batch = self.dPC.run_detection(self.neuronData['S'][i*batchSize:(i+1)*batchSize,:].toarray(),rerun=rerun,return_results=True,artificial=True)
        self.PC_fields = extend_dict(self.PC_fields,PC_fields_batch['fields']['parameter'].shape[0],PC_fields_batch)
        #return
        pickleData(PC_fields_batch,pathSave,'save')
  
  def assess_analysis(self,alpha=0.1,Z_thr=0):
    
    ## characterize results:
    ##% 1: pp = positive-positive: correctly detected PC
    ##% 2: ppf = positive-positive-false: correctly detected PC, but in wrong position
    ##% 3: pn = positive-negative: falsely detected PC
    ##% 4: np = negative-positive: not-detected PC
    ##% 5: nn = negative-negative: correctly not-detected PC
    
    nCells = self.PC_fields['fields']['parameter'].shape[0]
    print('adapt find_PCs, such that it can take these kind of values')
    #print(nCells)
    self.analysis = {'accuracy':np.zeros((nCells,5)).astype('bool'),
                'fr':np.zeros(nCells)}
    
    idx_MI = fdr_control(self.PC_fields['status']['MI_p_value'],alpha);
    idx_Bayes = self.PC_fields['status']['Bayes_factor'][:,0]-self.PC_fields['status']['Bayes_factor'][:,1] > Z_thr;
    
    PC_status = idx_MI & idx_Bayes
    self.PC_status = PC_status
    i = 0;
    #    accuracy = struct('pp',zeros(nCells,1),'ppf',zeros(nCells,1),'pn',zeros(nCells,1),'np',zeros(nCells,1),'nn',zeros(nCells,1));
    
    self.neuronPar = {'nModes':np.zeros(nCells)*np.NaN,
                      'theta':np.zeros(nCells)*np.NaN,
                      'sigma':np.zeros(nCells)*np.NaN,
                      'A':np.zeros(nCells)*np.NaN,
                      'A_0':np.zeros(nCells)*np.NaN,
                      'activation':np.zeros(nCells)*np.NaN}
    for n in range(nCells):
      for key in self.neuronData['fields'][n].keys():
        self.neuronPar[key][n] = self.neuronData['fields'][n][key]
    
    self.neuronPar['A_rate'] = self.neuronPar['A']/self.neuronPar['A_0']
    
    #print(neuronPar.keys())
    #print(neuronPar['nModes'].shape)
    
    theta_original = self.neuronPar['theta']
    theta_detected = self.PC_fields['fields']['parameter'][:,0,3,0]
    nbin = self.dPC.para['nbin']
    dTheta = (theta_detected - theta_original+nbin/2)%nbin-nbin/2
    
    for n in range(nCells):
      
      #if n/nCells>=perc*i:
        #print('>>> %3.0d%% done. time passed: %6.4g secs <<<'%(round(perc*i*100),time.time()-t_start))
        #i=i+1;
      
      #print('fields:')
      #print(PC_fields['fields']['parameter'][n,0,3,0])
      #print(neurons.field[n]['theta'])
      
      if PC_status[n] & (self.neuronPar['nModes'][n] > 0):#neuronData.PC_status[n]:
        if (self.PC_fields['fields']['nModes'][n] == self.neuronPar['nModes'][n]) & (abs(dTheta[n])<5):
        
          #print('save original & detected field position (including CI)')
        #((abs(PC_fields['fields']['parameter'][n,0,3,0] - neuronData['fields'][n]['theta'] + dPC.para['nbin']/2)%dPC.para['nbin']-dPC.para['nbin']/2) < 10):
          self.analysis['accuracy'][n,0] = True;
        else:
          self.analysis['accuracy'][n,1] = True;
        
      if PC_status[n] & ~(self.neuronPar['nModes'][n] > 0):#neuronData.PC_status[n]:
        self.analysis['accuracy'][n,2] = True;
  #        figure
  #        plot_ax_activity_artificial(gca,PC_fields(n),neuronData(n),bh,[])
  #        waitforbuttonpress
  #        close all
      elif ~PC_status[n] & (self.neuronPar['nModes'][n] > 0):#neuronData.PC_status[n]:
        self.analysis['accuracy'][n,3] = True;
  #        figure
  #        plot_ax_activity_artificial(gca,PC_fields(n),neuronData(n),bh,[])
  #        waitforbuttonpress
  #        close all
      elif ~PC_status[n] & ~(self.neuronPar['nModes'][n] > 0):#neuronData.PC_status[n]:
        self.analysis['accuracy'][n,4] = True;
  #        plt = false;
      
      #if plt_bool:
        #close all
        
        #plt.figure('position',[1000 500 1500 1200])
        #fields = plot_ax_activity_artificial(gca,PC_fields[n],neuronData(n),bh);
        
        #waitforbuttonpress
    
    #if nargin > 2:
      ### saving results
      #pathResults = pathcat(pathResults,sprintf('results_nC%d_s%d.mat',nCells,s))
      #save(pathResults,'analysis','PC_fields','neuronData','bh','-v7.3')
      
      #analyze_simulation(pathResults)
    
    
    #return neuronData, analysis, PC_fields
  def plot_analysis(self,key1,key2,key3,val3,nsteps=10):
    
    ratio = np.zeros((len(val3),nsteps,nsteps))*np.NaN
    N = np.zeros((len(val3),nsteps,nsteps))*np.NaN
    
    self.neuronPar['A_rate'] = self.neuronPar['A']/self.neuronPar['A_0']
    
    thresholds = {'A_0':np.linspace(0,10,nsteps+1),
                  'A_rate':np.linspace(0.25,4,nsteps+1),
                  'A':np.linspace(0,20,nsteps+1),
                  'activation':np.linspace(0,1,nsteps+1),
                  'sigma':np.linspace(1,8,nsteps+1),
                  'theta':np.linspace(0,100,nsteps+1)}
    #plt.figure(figsize=(8,3))
    #ax1 = plt.subplot(211)
    #ax2 = plt.subplot(212)
    #print(thresholds[key2])
    
    #ax1.plot([0,thresholds[key1][-1]],[0.2,0.2],'--',color=[0.8,0.8,0.8])
    #ax1.plot([0,thresholds[key1][-1]],[0.4,0.4],'--',color=[0.8,0.8,0.8])
    #ax1.plot([0,thresholds[key1][-1]],[0.6,0.6],'--',color=[0.8,0.8,0.8])
    #ax1.plot([0,thresholds[key1][-1]],[0.8,0.8],'--',color=[0.8,0.8,0.8])
    
    
    for k in range(len(val3)):
      for j in range(nsteps):
        for i in range(nsteps):
          idx = (self.neuronPar[key1] > thresholds[key1][i]) & (self.neuronPar[key1] <= thresholds[key1][i+1]) & (self.neuronPar[key2] > thresholds[key2][j]) & (self.neuronPar[key2] <= thresholds[key2][j+1]) & (self.neuronPar[key3] > val3[k])
          N[k,j,i] = idx.sum()
          ratio[k,j,i] = self.analysis['accuracy'][idx,0].sum()/N[k,j,i]
          
        #ax1.plot(thresholds[key1][:-1],ratio[k,j,:],label='%s > %5.3f'%(key2,thresholds[key2][j]))
        #ax2.plot(thresholds[key1][:-1],N[k,j,:],'o')
    #ax1.set_ylim([0,1])
    #ax1.set_xlim([thresholds[key1][0],thresholds[key1][-1]])  
    #ax2.set_xlim([thresholds[key1][0],thresholds[key1][-1]])  
    ##ax1.legend(loc='right')
    #plt.show(block=False)
    
    
    
    #f,axs = plt.subplots(2,2,figsize=(6,6))
    #idx = self.analysis['accuracy'][:,3]
    #axs[0][0].hist(self.neuronPar['A_0'][idx])
    #axs[0][0].set_xlabel('$A_0$')
    #axs[0][1].hist(self.neuronPar['A_rate'][idx])
    #axs[0][1].set_xlabel('$A_{rate}$')
    #axs[1][0].hist(self.neuronPar['activation'][idx])
    #axs[1][0].set_xlabel('activation')
    #axs[1][1].hist(self.neuronPar['sigma'][idx])
    #axs[1][1].set_xlabel('sigma')
    #plt.show(block=False)
    
    idxes = self.return_idxes(act_thr=0.3,A0_thr=0.5,Ar_thr=0.25)
    
    col_arr = ['g','r','k']
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(self.neuronPar['A_0'][idxes['false_neg']],self.neuronPar['A_rate'][idxes['false_neg']],self.neuronPar['activation'][idxes['false_neg']],s=2,color='r')
    ax.scatter(self.neuronPar['A_0'][idxes['true_pos']],self.neuronPar['A_rate'][idxes['true_pos']],self.neuronPar['activation'][idxes['true_pos']],s=2,color='b')
    ax.scatter(self.neuronPar['A_0'][idxes['false_theta']],self.neuronPar['A_rate'][idxes['false_theta']],self.neuronPar['activation'][idxes['false_theta']],s=2,color='k')
    ax.set_ylabel('$A_{0}$')
    ax.set_ylabel('$A_{rate}$')
    ax.set_zlabel('activation')
    plt.show(block=False)
    
    #plt.scatter(self.PC_fields['status']['MI_z_score'][idxes['false_neg']],self.PC_fields['status']['Bayes_factor'][idxes['false_neg'],0],s=2,color='r')
    #plt.scatter(self.PC_fields['status']['MI_z_score'][idxes['true_pos']],self.PC_fields['status']['Bayes_factor'][idxes['true_pos'],0],s=2,color='b')
    #plt.scatter(self.PC_fields['status']['MI_p_value'][idxes['false_neg']],self.neuronPar['A_rate'][idxes['false_neg']],s=2,color='r')
    #plt.scatter(self.PC_fields['status']['MI_p_value'][idxes['true_pos']],self.neuronPar['A_rate'][idxes['true_pos']],s=2,color='b')
    
    key_eval = ['MI_p_value','Bayes_factor']
    
    plt.figure()
    plt.subplot(121)
    plt.scatter(self.neuronPar['A_0'][idxes['thr_nPC'] & self.PC_status],self.PC_fields['status'][key_eval[0]][idxes['thr_nPC'] & self.PC_status],s=2,color='r')
    plt.scatter(self.neuronPar['A_0'][idxes['thr_nPC'] & ~self.PC_status],self.PC_fields['status'][key_eval[0]][idxes['thr_nPC'] & ~self.PC_status],s=2,color='g')
    plt.xlabel('A_0')
    plt.ylabel(key_eval[0])
    
    plt.subplot(122)
    plt.scatter(self.neuronPar['A_0'][idxes['thr_nPC'] & self.PC_status],self.PC_fields['status'][key_eval[1]][idxes['thr_nPC'] & self.PC_status,0],s=2,color='r')
    plt.scatter(self.neuronPar['A_0'][idxes['thr_nPC'] & ~self.PC_status],self.PC_fields['status'][key_eval[1]][idxes['thr_nPC'] & ~self.PC_status,0],s=2,color='g')
    plt.xlabel('A_0')
    plt.ylabel(key_eval[1])
    plt.show(block=False)
    #plt.xlabel('$A_0$')
    #plt.ylabel('Bayes_factor')
    #plt.title('non-place cells')
    
    id_eval = ['true_pos','false_neg','false_theta']
    plt.figure()
    
    for (key_y,j) in zip(key_eval,range(2)):
      for (key_x,i) in zip(['A_0','A','A_rate','activation'],range(4)):
        plt.subplot(2,4,j*4+i+1)
        for k in range(len(id_eval)):
          if len(self.PC_fields['status'][key_y].shape)>1:
            plt.scatter(self.neuronPar[key_x][idxes[id_eval[k]]],self.PC_fields['status'][key_y][idxes[id_eval[k]],0],s=2,color=col_arr[k])
          else:
            plt.scatter(self.neuronPar[key_x][idxes[id_eval[k]]],self.PC_fields['status'][key_y][idxes[id_eval[k]]],s=2,color=col_arr[k])
        plt.xlabel(key_x)
        plt.ylabel(key_y)
        plt.title('place cells')
    
    plt.show(block=False)
    
    
    plt.figure()
    plt.scatter(self.neuronPar['A_rate'][idxes['thr_PC']],self.PC_fields['fields']['parameter'][idxes['thr_PC'],0,1,0]/self.PC_fields['fields']['parameter'][idxes['thr_PC'],0,0,0])
    plt.plot([0,4],[0,4],'r')
    plt.ylim([0,10])
    plt.show(block=False)
    
    print('now, try inferring activation rate from errorbars of A')
    print('use average amplitude to infer activity baseline?!')
    
    plt.figure()
    plt.scatter(self.neuronPar['activation'][idxes['thr_PC']],self.PC_fields['fields']['parameter'][idxes['thr_PC'],0,1,1])#/self.PC_fields['fields']['parameter'][idxes['thr_PC'],0,0,0])
    plt.plot([0,1],[0,1],'r')
    plt.ylim([0,10])
    plt.show(block=False)
    
    
    X, Y = np.meshgrid(thresholds[key1][:-1], thresholds[key2][:-1])
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for k in range(len(val3)):
      ax.plot_surface(X,Y,ratio[k,:,:],cmap=cm.jet)
    ax.set_xlabel(key1)
    ax.set_ylabel(key2)
    ax.set_zlabel('performance')
    plt.show(block=False)
    
  def RoC(self,zrange=[-10,10],nsteps=21):
    
    tp = np.zeros(nsteps)
    fp = np.zeros(nsteps)
    ftheta = np.zeros(nsteps)
    
    idxes = self.return_idxes(act_thr=0.7,A0_thr=0.5,Ar_thr=0.25)
    
    N_PC = (idxes['thr_PC'] & (self.neuronPar['nModes']>0)).sum()
    N_nPC = (idxes['thr_nPC'] & (self.neuronPar['nModes']==0)).sum()
    
    print('place cells: %d'%N_PC)
    print('non-place cells: %d'%N_nPC)
    
    z_arr = np.linspace(zrange[0],zrange[1],nsteps)
    for i in range(nsteps):
      self.assess_analysis(alpha=1,Z_thr=z_arr[i])
      idxes = self.return_idxes(act_thr=0.7,A0_thr=0.5,Ar_thr=0.25)
      
      tp[i] = self.analysis['accuracy'][idxes['true_pos'],0].sum()/N_PC
      fp[i] = self.analysis['accuracy'][idxes['false_pos'],2].sum()/N_nPC
      ftheta[i] = self.analysis['accuracy'][idxes['false_theta'],1].sum()/N_PC
    
    Z_thr = 0
    Z_idx = np.argmin(abs(z_arr-Z_thr))
    plt.figure()
    plt.subplot(122)
    plt.plot([0,1],[0,1],'--',color=[0.6,0.6,0.6])
    plt.plot(fp,tp,'k')
    plt.plot(fp[Z_idx],tp[Z_idx],'ro')
    plt.xlim([0,1])
    plt.ylim([0,1])
    
    plt.subplot(121)
    plt.plot(z_arr,tp,'g')
    plt.plot(z_arr,fp,'r')
    plt.plot(z_arr,ftheta,'b')
    
    plt.plot(z_arr[Z_idx],tp[Z_idx],'ko')
    plt.plot(z_arr[Z_idx],fp[Z_idx],'ko')
    plt.plot(z_arr[Z_idx],ftheta[Z_idx],'ko')
    plt.ylim([0,1])
    
    plt.show(block=False)
    
    
  
  def return_idxes(self,act_thr=0.7,A0_thr=1,Ar_thr=0.25):
    idx_thr_PC = (self.neuronPar['activation']>act_thr) & (self.neuronPar['A_0']>A0_thr) & (self.neuronPar['A_rate']>Ar_thr)
    idx_thr_nPC = (self.neuronPar['A_0']>A0_thr)
    idxes = {'thr_PC':idx_thr_PC,
              'thr_nPC':idx_thr_nPC,
              'true_pos':self.analysis['accuracy'][:,0] & idx_thr_PC,
              'false_theta':self.analysis['accuracy'][:,1] & idx_thr_PC,
              'false_pos':self.analysis['accuracy'][:,2] & idx_thr_nPC,
              'false_neg':self.analysis['accuracy'][:,3] & idx_thr_PC,
              'true_neg':self.analysis['accuracy'][:,4] & idx_thr_nPC}
    return idxes
#nP = 4
#nCells = 100
#AD = artificialData("/media/wollex/Analyze_AS1/linstop","762",10,nP,redo=True)
#AD.generate_activity(nCells)
#AD.analyze(batchSize=nCells/2)

#PCs, neurons = generate_artificial_data("/media/wollex/Analyze_AS3/Data","879",7,nCells,nP)
#PCs, neurons = generate_artificial_data("/media/wollex/Analyze_AS3/Data","879",10,nCells,nP)


## do analysis: which PCs are best detected (width, amplitude, firingrate, ...), number of fields, etc
## enable multiple fields
## baseline shouldn't be monotonous, rather random distributed?!

  