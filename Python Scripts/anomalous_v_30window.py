import netCDF4 as nc
import numpy as np
import os
import math
import matplotlib.pyplot as plt

#def standardization(data):
#    mu = np.mean(data, axis=0)
#    sigma = np.std(data, axis=0)
#    print(mu, sigma)
#    return (data - mu) / sigma

# yy = smooth(y) smooths the data in the column vector y ..
# The first few elements of yy are given by
# yy(1) = y(1)
# yy(2) = (y(1) + y(2) + y(3))/3
# yy(3) = (y(1) + y(2) + y(3) + y(4) + y(5))/5
# yy(4) = (y(2) + y(3) + y(4) + y(5) + y(6))/5
# ...

def smooth(a, WSZ):
    # a:原始数据，NumPy 1-D array containing the data to be smoothed
    # 必须是1-D的，如果不是，请使用 np.ravel()或者np.squeeze()转化 
    # WSZ: smoothing window size needs, which must be odd number,
    # as in the original MATLAB implementation
    out0 = np.convolve(a, np.ones(WSZ,dtype=int),'valid')/WSZ
    r = np.arange(1, WSZ-1, 2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((  start , out0, stop  ))

def insertzeros(t, x, zero=0):
    ta = []
    positive = (x-zero) > 0
    ti = np.where(np.bitwise_xor(positive[1:], positive[:-1]))[0]
    for i in ti:
        y_ = np.sort(x[i:i+2])
        z_ = t[i:i+2][np.argsort(x[i:i+2])]
        t_ = np.interp(zero, y_, z_)
        ta.append( t_ )
    tnew = np.append( t, np.array(ta) )
    xnew = np.append( x, np.ones(len(ta))*zero )
    xnew = xnew[tnew.argsort()]
    tnew = np.sort(tnew)
    return tnew, xnew

root_Path = r'I:\ERA5\vwnd_1979_2019'
fileLists = os.listdir(root_Path) #need to change!!!!!!!!!!!!!!!!!!!!!!!!
ds = nc.Dataset(r'I:\ERA5\vwnd_1979_2019\vwnd.1979.nc')
var = 'v'
nlon = 11
nlat = 10
nyear = 41
nday = 365

v = np.empty((nyear, nday, nlat, nlon), dtype=float)
for i, file in enumerate(fileLists):
    #var_info = ds.variables[var]
    file_name = os.path.join(root_Path, file)
    ds = nc.Dataset(file_name)
    var_data = ds[var][:, 58:48:-1, 110:121]
    if len(var_data[:, 0, 0]) == 1464: #whether it's a leap year?
        var_data = np.delete(var_data, (236, 237, 238, 239), axis = 0)
    var_data = np.array(var_data).reshape(365, 4, nlat, nlon)
    var_data_day = np.mean(var_data, 1) #4th average
    v[i, :, :, :] = var_data_day
    print(file)
lat = ds.variables["latitude"][58:48:-1]
#var_data_cal_annual = np.mean(var_data_cal, ) #getting the daily average

### calculate the 31 years average and get the anomaly
v_avg = np.mean(v, axis=0)

## running mean
run = np.empty((nday, nlat, nlon))
for i in range(0, nlat):
    for j in range(0, nlon):
        run[:, i, j] = smooth(v_avg[:, i, j], 31)

## the anomaly
v_ano = np.empty((nyear, nday, nlat, nlon))
for i in range(0, nyear):
    v_ano[i, :, :, :] = v[i, :, :, :] - run
## then doing 7 days running mean for anomaly
v_ano = v_ano.reshape((nyear*nday, nlat, nlon))
v_ano_run = np.empty(v_ano.shape)
for i in range(0, nlat):
    for j in range(0, nlon):
        v_ano_run[:, i, j] = smooth(v_ano[:, i, j], 7)
del v_ano, run, v_avg
### the area average for the North China
weights = [math.cos(i / 180. * math.pi) for i in lat]
v_area = []
for k in range(0, nday*nyear):
    avg = 0.
    for j in range(0, nlat):
        for i in range(0, nlon):
            avg = avg + v_ano_run[k, j, i]* weights[j] / float(nlat) / float(nlon)
    v_area = np.append(v_area, avg)

### fetch the winter data from November of the first year to March of the second year
v_cal = []
for i in range(0, nyear-1):
    v_cal = np.append(v_cal, v_area[i * 365 + 304 : (i + 1) * 365 + 89 + 1])
print(v_cal.shape)
num_day = len(v_cal)
  
critical_std1 = -0.
### start to find the window that anomalous v beyond 30 days 
lasting = []
begin = []
i = 15
while i < num_day:
    if (int((i + 15) / 151) - int((i - 15) / 151) != 0):  
        i = i + 1
        continue         
    
    tempor = v_cal[i - 15 : i + 15]
    ano_day = len(v_cal[np.where(tempor <= critical_std1)]) # find the index of the anomalous v 
    ano_index = []
    while ano_day >= 0.75 * 31:
        ano_index = np.append(ano_index, i)
        i = i + 1
        tempor = v_cal[i - 15 : i + 15]
        ano_day = len(v_cal[np.where(tempor <= critical_std1)]) # find the index of the anomalous v 
    if len(ano_index) != 0:
        begin = np.append(begin, ano_index[0] - 15)
        lasting = np.append(lasting, ano_index[-1] - ano_index[0] + 30 + 1)
        i = i + 14 + 15
        continue
    i = i + 1
   
np.savetxt(r'I:\ERA5\anomalous_v_30window_try.txt', \
    np.column_stack((lasting, begin)), fmt='%d', newline = '\n', delimiter = ' ')

total_event = len(begin) #the numbers of anomalous v component of winds events
max_ind = np.full(total_event, 0)
for i in range(0, total_event):
    tempor = v_cal[int(begin[i]) : int(begin[i] + lasting[i] - 1)]
    max_ind[i] = np.where(tempor == min(tempor))[0][0] + begin[i]

#fig, axs = plt.subplots(2, 1)
fig = plt.figure(figsize=(9, 8))
t = np.arange(-30, 31, 1)
lower = 0.
for i in range(total_event):
    ax = fig.add_subplot(6, 7, i+1)
    tempor = v_cal[-30 + max_ind[i] : 31 + max_ind[i]]
    #ax.xaxis.set_major_locator(range(-7, 10, 1))
    t1,x1 = insertzeros(t, tempor, zero=lower)

    xm = np.copy(x1)
    xm[x1 < lower] = np.nan        
    ax.plot(t1, xm, color="#fd5f00", linewidth=1)

    xl = np.copy(x1)
    xl[(x1 > lower)] = np.nan        
    ax.plot(t1, xl, color="#005792", linewidth=1)
    ymax = 4
    ymin = -4
    ax.text(-25, ymax+1, str(int(max_ind[i] - begin[i] + 1)), fontsize = 'x-small',
            ha='center', va='top', weight='bold', color='black')
    ax.text(25, ymax+1, str(int(lasting[i])), fontsize = 'x-small',
            ha='center', va='top', weight='bold', color='black')

    ax.tick_params(axis='both', which='major', pad=1.5)
    ax.set_yticks(range(ymin, ymax, 1))
    ax.set_yticklabels(range(ymin, ymax, 1), fontsize = 3.5)
    ax.set_xticks(range(-25, 35, 5))
    ax.set_xticklabels(range(-25, 30, 5), fontsize = 3.5)

    ax.set_xlim(-30, 30)
    ax.set_ylim(ymin, ymax)
    #ax.set_xlabel('time')
    #ax.set_ylabel('v wind')
    ax.grid(True)
plt.subplots_adjust(wspace = 0.2, hspace = 0.3)
fig.savefig("E:\\yoghuit\\Desktop\\anomalous_v.jpeg",dpi = 600)

plt.show()