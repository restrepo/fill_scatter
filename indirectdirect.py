
#How to compute sigmav_bb, sigmav_w+w-, LUX y XENON1T(2017) constraints

''' Experimental annihilation cross section, lower limit taken
from b bbar and W+W- as a function of the dark matter mass, dawarf galaxies 
Fermi limit 2014 http://arxiv.org/pdf/1503.02641.pdf
'''

import numpy as np
def sigmav_bb_fit(x):

    c = 2.15034522e-32
    b = 1.66066271e-28
    a = 2.59843102e-27
    
    g = a + b*x + c*x**2 
    return g  

def sigmav_ww_fit(x):

    c = 6.16045106e-32
    b = 1.88303853e-28
    a = 1.08024182e-26
    
    g = a + b*x + c*x**2 
    return g  

def LUXconstraint_official(x):
    '''Evaluate LUX limits for DD [arXiv: astro-ph.CO/1310.8214 ]'''
    m = np.array([7.89119649, 10., 12., 17., 30., 200., 450., 1000., 2000., 3500., 5000.])
    A = np.array([0.00000000e+00, 8.48911168e-41, 1.95004609e-42, 3.75121862e-44, 3.28089421e-45, 5.80846495e-46, 1.39568112e-45, 3.05245738e-45, 5.63013735e-45, 1.07867168e-44, 1.93910180e-44])
    B = np.array([0.00000000e+00, -9.02408269e-01, -5.20772861e-01, -1.90971186e-01, -5.03032503e-02, 7.66657692e-03, 3.11852240e-03, 1.34049663e-03, 7.28957975e-04, 3.94487035e-04, 2.23523977e-04])
    imax=10
    for i in range(1,imax+1):
        if x>=m[i-1] and x<=m[i]:
            limit = A[i]*np.exp(B[i]*x)
        if x<m[0] or x>m[imax]:
            limit = 0.0	
            print 'ERROR: Out of range'
    return limit*1E36

def sigmaSI_LUX_f(x):
    import numpy as np
    x=np.asarray(x)
    if not x.shape: 
        x=[x] 
    
    return map(LUXconstraint_official,x)

def XENON1Tconstraint_official(x):
    '''Evaluate XENON1T limits for DD [arXiv: 1206.6288v1]'''
    m = np.array([5.13300657, 6., 7., 8., 10., 12., 15., 22., 28., 40., 70., 170., 250., 400., 900., 2000.,3100 ])
    A = np.array([0.00000000e+00, 9.12202809e-31, 1.05982015e-34, 6.87776035e-38, 7.11912593e-40, 9.09058495e-42, 4.75346716e-43, 8.00154523e-45, 3.48707854e-46, 1.06581517e-46, 2.23043755e-47, 1.24518124e-47, 1.53905355e-47, 2.05705600e-47, 3.78566641e-47, 7.31147326e-47, 1.27722554e-46 ])
    B = np.array([0.00000000e+00, -4.27599323e+00, -2.76279703e+00, -1.71395963e+00, -1.16210971e+00, -7.44082481e-01, -5.02777951e-01, -2.35307322e-01, -8.70057359e-02, -4.52449357e-02, -3.62262032e-03, 5.56327336e-03, 4.32584727e-03, 3.10054027e-03, 1.60037074e-03, 8.62949144e-04, 6.50844801e-04 ])
    imax=15
    for i in range(1,imax+1):
        if x>=m[i-1] and x<=m[i]:
            limit = A[i]*np.exp(B[i]*x)
        if x<m[0] or x>m[imax]:
            limit = 0.0	
            print 'ERROR: Out of range'
    return limit*1E36

##Sigma en pb (1cm^{2}=10^{-36})
def sigmaSI_XENON1T_f(x):
    import numpy as np
    x=np.asarray(x)
    if not x.shape: 
        x=[x] 
    
    return map(XENON1Tconstraint_official,x)