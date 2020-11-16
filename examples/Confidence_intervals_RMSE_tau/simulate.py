import numpy as np
import scipy.stats
import bootstrapped.bootstrap as bs

def RMSE_function(values):
  tmp = np.sqrt(np.mean(values*values,axis=1))
  return tmp


def tau_function(a):
  s = a.shape
  tmp = np.array([])
  if len(s) == 2:
    tmp = scipy.stats.kendalltau(a[:,0],a[:,1])
  else:
    for i in range(s[0]):
      x = scipy.stats.kendalltau(a[i,:,0],a[i,:,1])[0]
      if x == np.nan:
        tmp = np.append(tmp,0.0)
      else:
        tmp = np.append(tmp,x)

  return tmp


def eval_stats(exp,pred,N):
  error = np.abs(exp-pred)
  rmse_bootstrap_dist = bs.bootstrap(np.reshape(error,(N,-1)), stat_func=RMSE_function,num_iterations=1000, alpha =0.05,is_pivotal=True,return_distribution=True)
  rmse = np.sqrt(np.mean(error**2))
  tau = scipy.stats.kendalltau(exp,pred)[0]
  tau_bootstrap_dist = bs.bootstrap(np.reshape(np.array(list(zip(exp,pred))),(N,-1)), stat_func=tau_function,num_iterations=1000, alpha =0.05,is_pivotal=True,return_distribution=True)
  tmp1 = [rmse,np.percentile(rmse_bootstrap_dist,2.5),np.percentile(rmse_bootstrap_dist,97.5),tau,np.percentile(tau_bootstrap_dist,2.5),np.percentile(tau_bootstrap_dist,97.5)]
  return tmp1


#generate toy data in interval [-12:-5]
all_results = []
NN = [10,25,35,50,75,100,200]
for N in NN:
#  simulate experimental data as uniform distribution
  data = [ np.random.uniform(low=-12,high=-5,size=N) for x in range(1000)]
  tmp = []
  for i in data:
    #simulated predicted data with gaussian error and calculate statistics for datasets
    tmpresults = []
    for sigma in [0.5,1.0,1.5,2.0]:
      i_pred = i + np.random.normal(0,sigma,N)
      results = eval_stats(i,i_pred,N)
      tmpresults.append(results)

    tmp.append(tmpresults)
    
  
  tmp = np.array(tmp)
  tmp = np.reshape(tmp,(1000,24))
  all_results.append(np.mean(tmp,axis=0))

all_results = np.array(all_results)
all_results = np.reshape(all_results,(-1,24))
import pandas as pd
df = pd.DataFrame(all_results,columns=['mean_rmse_0.5','low_rmse_0.5','high_rmse_0.5','mean_tau_0.5','low_tau_0.5','high_tau_0.5','mean_rmse_1.0','low_rmse_1.0','high_rmse_1.0','mean_tau_1.0','low_tau_1.0','high_tau_1.0','mean_rmse_1.5','low_rmse_1.5','high_rmse_1.5','mean_tau_1.5','low_tau_1.5','high_tau_1.5','mean_rmse_2.0','low_rmse_2.0','high_rmse_2.0','mean_tau_2.0','low_tau_2.0','high_tau_2.0'])
df['N']= NN
df.to_csv('results_N_CI_error.csv')

