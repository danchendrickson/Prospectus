import numpy as np
#import matplotlib.pyplot as plt

import mpu6050 as mp
import time
import zipfile
import pandas as pd
import os
path = '/'


sensor = mp.mpu6050(0x68)

start = time.time()

j = 1

GB = 2

while GB > 1:
    Results = []
    FileStart = time.time()
    for i in range(60000):
        gx, gy, gz = sensor.get_accel_data(True)
        #tx, ty, tz = sensor.get_gyro_data()
        Results.append([float(int((time.time()-FileStart) * 1000.0))/1000.0, gx, gy, gz])
        time.sleep(0.0001)

    Results = np.matrix(Results)

    NextFileName = '/home/pi/AccelData/60kPoints-'+time.strftime('%y%m%d-%H%M')
    
    df = pd.DataFrame(data=Results.astype('float'))
    df.to_csv(NextFileName+'.csv', sep=',', header=False, float_format='%.5f')
    #with open(NextFileName+'.csv','wb') as f:
    #    for line in Results:
    #        np.savetxt(f, line, fmt='%.5f')
    #Results.tofile('/home/pi/Code/Prospectus/Data/Accel/60kPoints'+str(j)+'.txt', sep=',')
    print(j, (time.time()-start) / 60.0)
    j+=1

    zip_file = zipfile.ZipFile(NextFileName+'a.zip', 'w')
    zip_file.write(NextFileName+'.csv', compress_type=zipfile.ZIP_DEFLATED)
    zip_file.close()
    
    os.remove(NextFileName+'.csv')
    
    st = os.statvfs(path)
    bytes_avail = (st.f_bavail * st.f_frsize)
    GB = bytes_avail / 1024 / 1024 / 1024