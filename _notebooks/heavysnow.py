import numpy as np

class GraphSignal: 
    def __init__(self,V,W,f):
        self.f=np.array(f) 
        self.V=V 
        self.W=np.array(W)
        self.n=len(self.f)
        self.degree=np.sum(np.array(self.W),0)
        self.initdist=self.degree/np.sum(self.degree) 
        
class HeavySnowTransform:  
    def __init__(self,graphsignal):
        self.n=graphsignal.n
        self.f=graphsignal.f
        self.V=graphsignal.V
        self.nodeindex=np.arange(0,self.n,dtype='int')
        self.graphweight=graphsignal.W
        self.initdist=graphsignal.initdist
        self._initialize()
    
    def _initialize(self):
        self.tau=0
        self.b=0
        self.trajectory=[np.random.choice(self.n,1,p=self.initdist).item()]
        self.flowcount=[0]
        self.snowygrounds=None 
        self.eucliddistance=(np.array(self.f)[:,np.newaxis]-np.array(self.f)[np.newaxis,:])**2
        self.snowdistance=np.zeros((self.n,self.n))
        self.euclidweight=np.zeros((self.n,self.n))
        self.snowweight=np.zeros((self.n,self.n))
        
    def _snowonce(self,ell,maxflow):
        b=self.b
        flowcount=self.flowcount[-1] 
        snowyground=self.snowygrounds[...,ell-1].copy()
        currentnode=self.trajectory[-1] 
        neighbor=self.nodeindex[(self.graphweight>0)[currentnode]]
        transitionprob=None
        nextnode=None 
        downstream=None
        
        if flowcount > maxflow: # reset 
            nextnode=np.random.choice(self.n,1,p=self.initdist).item()
            snowyground[nextnode]=snowyground[nextnode]+b
            flowcount=0
        elif sum(neighbor)==0: # empty neighborhood 
            nextnode=np.random.choice(self.n,1,p=self.initdist).item()
            snowyground[nextnode]=snowyground[nextnode]+b
            flowcount=0
        else:
            _transitionprob=self.graphweight[currentnode]/sum(self.graphweight[currentnode])
            nextnode=np.random.choice(self.n,1,p=_transitionprob).item()
            downstream=neighbor[np.where(snowyground[currentnode]>=snowyground[neighbor])]
            if nextnode in set(downstream): # flow 
                snowyground[nextnode]=snowyground[nextnode]+b
                flowcount=flowcount+1
            else: # block 
                nextnode=np.random.choice(self.n,1,p=self.initdist).item()
                snowyground[currentnode]=snowyground[currentnode]+b
                flowcount=0
        self.snowygrounds[...,ell]=snowyground
        self.flowcount=self.flowcount+[flowcount]
        self.trajectory=self.trajectory+[nextnode]
        
    def _getdegreematrix(self,matrix):
        return np.diag(np.sum(matrix,1))
    
    def _normalize(self,matrix): 
        return np.sqrt(self._getdegreematrix(matrix))@matrix@np.sqrt(self._getdegreematrix(matrix))
    
    def _updateeuclidweight(self):
        self.euclidweight=np.exp(-self.eucliddistance/(self.b**2))-np.eye(self.n)
        
    def _updatesnowdistance(self):
        self.snowdistance=np.sum((self.snowygrounds[:,np.newaxis,:]-self.snowygrounds[np.newaxis,:,:])**2,axis=-1)
             
    def _updatesnowweight(self):
        self.snowweight=np.exp(-self.snowdistance/(self.tau*self.b**2))-np.eye(self.n)
    
    def snow(self,tau,b=1,maxflow=99999):
        self._initialize()
        self.b=b
        self.tau=tau
        self.snowygrounds=np.repeat(self.f,self.tau+1).reshape(self.n,self.tau+1)
        print('HST (tau= %s, b=%s)' % (self.tau,self.b))
        for ell in np.arange(1,self.tau+1,dtype='int'): 
            print('\r'+str(ell)+'/'+str(self.tau),sep='',end='')
            self._snowonce(ell,maxflow)
        print('\n'+'HST completed and all history is recorded.')
        self._updatesnowdistance()
        self._updatesnowweight()
        self._updateeuclidweight()                
        
# class HeavySnowTransform_old:  #hst(f,W,V,τ=τ,b=0.03,γ=1)
#     def __init__(self,graphsignal):
#         self.n=graphsignal.n
#         self.f=graphsignal.f
#         self.V=graphsignal.V
#         self.graphweight=graphsignal.W
#         #self.degree=graphsignal.degree
#         self.initdist=graphsignal.initdist
#         self._initialize()
    
#     def _initialize(self):
#         self.tau=0
#         self.b=0
#         self.trajectory=[0]
#         self.flowcount=[0]
#         self.snowygrounds=pd.DataFrame(data={'h0':self.f})            
#         self.snowygrounds_f=pd.DataFrame(data={'h0':self.f})            
#         self.snowygrounds_sumxi=pd.DataFrame(data={'h0':self.f*0})            
#         self.eucliddistance=(np.array(self.f)[:,np.newaxis]-np.array(self.f)[np.newaxis,:])**2
#         self.snowdistance=np.zeros((self.n,self.n))
# #         self.snowdistance1=np.zeros((self.n,self.n))
# #         self.snowdistance2=np.zeros((self.n,self.n))
# #         self.snowdistance3=np.zeros((self.n,self.n))
#         self.euclidweight=np.zeros((self.n,self.n))
#         self.snowweight=np.zeros((self.n,self.n))
        
#     def _snowonce(self,ell,b,maxflow):
#         W=self.graphweight>np.asmatrix(np.random.random((self.n,self.n)))
#         hnext=np.array(self.snowygrounds.iloc[:,-1])
#         currentnode=self.trajectory[-1]
#         flowcount=self.flowcount[-1]
#         # 1. h(u) <- h(u)+b : update current node when flowcound > 0 
#         hnext[currentnode]=hnext[currentnode]+b
#         # 2. check that: are there any nodes to which snow can flow from u. 
#         neighbor=np.where(W[currentnode,:]>0)[1] ## Nu is np.array
#         if len(neighbor)==0: downstream=np.array([])
#         else: downstream=neighbor[list((np.where(hnext[neighbor]<=hnext[currentnode]))[0])]
#         # 3. check flowable 
#         flowable = len(downstream)>0 and flowcount < maxflow 
#         # 4. determine flow or block
#         if flowable==0: # block!
#             node0=np.random.choice(self.n,1,p=self.initdist).item()
#             neighbor=np.where(W[node0,:]>0)[1] 
#             downstream=neighbor[list((np.where(hnext[neighbor]<=hnext[node0]))[0])]
#             if len(downstream)==0: 
#                 nextnode=node0
#             else: 
#                 nextnode=np.random.choice(list(downstream),1).item()
#             flowcount=0
#         else: #flow
#             nextnode=np.random.choice(list(downstream),1).item()
#             flowcount=flowcount+1
#         self.snowygrounds['h%s' % ell]=hnext # pd.concat 
#         self.flowcount=self.flowcount+[flowcount] #pd.concat 
#         self.trajectory=self.trajectory+[nextnode] #pd.concat 
        
#         self.snowygrounds_f['h%s' % ell]=self.f
#         self.snowygrounds_sumxi['h%s' % ell]=hnext-self.f

#     def _getdegreematrix(self,matrix):
#         return np.diag(np.sum(matrix,1))
    
#     def _normalize(self,matrix): 
#         return np.sqrt(self._getdegreematrix(matrix))@matrix@np.sqrt(self._getdegreematrix(matrix))
    
#     def _updateeuclidweight(self):
#         self.euclidweight=np.exp(-self.eucliddistance/(self.b**2))-np.eye(self.n)
        
#     def _updatesnowdistance(self):
#         self.snowdistance=np.sum((np.array(self.snowygrounds)[:,np.newaxis,:]-np.array(self.snowygrounds)[np.newaxis,:,:])**2,axis=-1)
        
# #     def _updatesnowdistance1(self):
# #         self.snowdistance1=np.sum((np.array(self.snowygrounds_f)[:,np.newaxis,:]-np.array(self.snowygrounds_f)[np.newaxis,:,:])**2,axis=-1)
        
# #     def _updatesnowdistance2(self):
# #         self.snowdistance2=np.sum(2*(np.array(self.snowygrounds_f)[:,np.newaxis,:]-np.array(self.snowygrounds_f)[np.newaxis,:,:])*(np.array(self.snowygrounds_sumxi)[:,np.newaxis,:]-np.array(self.snowygrounds_sumxi)[np.newaxis,:,:])   ,axis=-1)

# #     def _updatesnowdistance3(self):
# #         self.snowdistance3=np.sum((np.array(self.snowygrounds_sumxi)[:,np.newaxis,:]-np.array(self.snowygrounds_sumxi)[np.newaxis,:,:])**2,axis=-1)
        
#     def _updatesnowweight(self):
#          self.snowweight=np.exp(-self.snowdistance/(self.tau*self.b**2))-np.eye(self.n)
#     def snow(self,tau,b=1,maxflow=100):
#         self._initialize()
#         self.b=b
#         self.tau=tau
#         print('HST (tau= %s, b=%s)' % (tau,b))
#         for ell in np.arange(1,tau+1): 
#             print('\r'+str(ell)+'/'+str(tau),sep='',end='')
#             self._snowonce(ell,b,maxflow)
#         print('\n'+'HST completed and all history is recorded.')
#         self._updatesnowdistance()
#         self._updatesnowdistance1()
# #        self._updatesnowdistance2()
# #        self._updatesnowdistance3()
# #        self._updatesnowweight()
#         self._updateeuclidweight()