#!/usr/bin/env python
# coding: utf-8

# In[31]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# In[32]:


isotopes = np.array(["Sb122" , "Sb124","Al28","Cu64","Cu66"])
t12_day  = np.array([ 2.7238 , 60.20, 0.0015625, 0.529167, 0.003556])


# In[33]:


mcnp_outp = np.array([8.28920E-06,4.89742E-06,5.95396E-07,2.61495E-07,2.75076E-07])
mcnp_unc  = np.array([ 1.26, 2.31, .42, 0.6, 0.48])


# In[34]:


## Simulation settings
vol_cc = np.array([0.14945,0.14945,14.37755,14.37755,14.37755])
nps    = 1e7


# In[35]:


data = {'isotope':isotopes, 't12_day':t12_day, 
        'mcnp_outp':mcnp_outp, 'vol':vol_cc}
 
# Create DataFrame
df = pd.DataFrame.from_dict(data)
df["t12_hour"] = 24. * df["t12_day"]
 
# Print the output.
print(df)


# In[36]:


df["prod-per-n"] = nps * df["vol"] * df["mcnp_outp"]


# In[37]:


dd_rate = 1e8 ## max rate of n/sec


# In[38]:


## In one second, DD creates 10^8 neutrons
sims_per_sec = dd_rate/nps
print(sims_per_sec)


# In[39]:


exposure_hrs = 90.0
df["max-n-activated"] = (exposure_hrs*3600)*df["prod-per-n"]*sims_per_sec


# In[ ]:





# In[40]:


df["after-act-prof"] = df["max-n-activated"] * (1-np.exp(-np.log(2) *(exposure_hrs/24.) /df["t12_day"]))*( (1/np.log(2))*(df["t12_day"]/(exposure_hrs/24.)))


# In[41]:


cooldown_days = 20
df["after-decay-time"] = df["after-act-prof"] * np.exp(-np.log(2)*cooldown_days/df["t12_day"])


# In[42]:


df['act-Bq'] = df["after-decay-time"]*np.log(2)/(24*3600*df["t12_day"])
df['act-kBq']=df['act-Bq']/1000


# In[43]:


print(df)


# In[44]:


def line(x,m,b):
    return m*x+b
from scipy.optimize import curve_fit
x_data=np.array([range(1,100)])
def calcActivity(exposure_hrs):
    calcAct = (x_data*3600)*7.319194*sims_per_sec
    calcAct = calcAct * (1-np.exp(-np.log(2) *(exposure_hrs/24) /60.2))*( (1/np.log(2))*(60.2/(exposure_hrs/24)))
    calcAct = calcAct * np.exp(-np.log(2)*20/60.2)
    calcAct = (calcAct *np.log(2)/(24*3600*60.2))/1000
    return calcAct
y_data=calcActivity(x_data)
x_data.shape = (99,)
y_data.shape = (99,)
plt.figure(2)
plt.scatter(x_data,y_data, label='Activation Data')
plt.xlabel('Exposure time in hours')
plt.ylabel('Final activity of Sb124 in kBq')


param_guess = np.array([ 1e-7 , 1e-7])
gss_y_vals = line(x_data, param_guess[0], param_guess[1],)

# -- Do the fit, plot it
popt, pcov = curve_fit(line, x_data, y_data, p0=param_guess)
perr = np.sqrt(np.diag(pcov))
fit_y_vals = line(x_data, popt[0], popt[1],)
plt.plot(x_data, fit_y_vals, color='r', label="Best Line fit")
print(popt)
print(perr)

# -- Turn on a plot key
plt.legend(loc='upper left')


# In[45]:


def curve(x,a,c):
    return a*x**-2 + c
from scipy.optimize import curve_fit
x_data=np.array([26.5,38.5,50.5,62.5,74.5])
y_data=np.array([1.323492e-03,2.141463e-04,7.639059e-05,1.344793e-05,8.638147e-06])
x_data.shape = (5,)
y_data.shape = (5,)
plt.figure(2)
plt.scatter(x_data,y_data, label='Activation Data')
plt.xlabel('Antimony distance from neutron source in centimeters')
plt.ylabel('Final activity of Sb124 in kBq')


param_guess = np.array([ 50, 5])
gss_y_vals = curve(x_data, param_guess[0], param_guess[1])

# -- Do the fit, plot it
popt, pcov = curve_fit(curve, x_data, y_data, p0=param_guess)
perr = np.sqrt(np.diag(pcov))
fit_y_vals = curve(x_data, popt[0], popt[1])
plt.plot(x_data, fit_y_vals, color='r', label="Best Line fit")
print(popt)
print(perr)
# -- Turn on a plot key
plt.legend(loc='upper right')


# In[ ]:




