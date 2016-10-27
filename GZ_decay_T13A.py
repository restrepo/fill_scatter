
# It takes a pandas.core.frame.DataFrame something like: 
#y=xd[['Mchi1','Mchi2','Mchi3','M_DM','N21','N22','N23','N31','N32','N33']]

import numpy as np
import pandas as pd

def z_decay(y):
    
    Z_Decay=[]
    y=np.asarray(y)
    
    MZ = 91.1876
    gw = 0.6483971671950267
    cw = 0.8753861242778386   
    
    for i in range(0,len(y)):
        #M_DM=y[i,3]
        if MZ > 2*y[i,3]:
            
            if y[i,3] == np.abs(y[i,0]):
                
                gammaZ = MZ/(6.*np.pi)*(1.-4*(y[i,3]/MZ)**2)**(3./2.)*(gw/(4.*cw)*(y[i,7]**2-y[i,4]**2))**2  
                
            if y[i,3] == np.abs(y[i,1]):
                
                gammaZ = MZ/(6.*np.pi)*(1.-4*(y[i,3]/MZ)**2)**(3./2.)*(gw/(4.*cw)*(y[i,8]**2-y[i,5]**2))**2   
                
            if y[i,3] == np.abs(y[i,2]):
                
                gammaZ = MZ/(6.*np.pi)*(1.-4*(y[i,3]/MZ)**2)**(3./2.)*(gw/(4.*cw)*(y[i,9]**2-y[i,6]**2))**2                   
        else:
            gammaZ = 0.
        
        Z_Decay.append([gammaZ])
    
    zd=pd.DataFrame(Z_Decay,columns=['Z_Decay'])
    
    return  zd