import numpy as np
import matplotlib.pyplot as plt

import mpu6050 as mp
import time
import zipfile
import pandas as pd

sensor = mp.mpu6050(0x68)

start = time.time()

j = 1

while ((time.time() - start) / 60.0 /60.0) < 0.125:
    Results = []
    for i in range(60000):
        gx, gy, gz = sensor.get_accel_data(True)
        #tx, ty, tz = sensor.get_gyro_data()
        Results.append([float(int((time.time()-start) * 1000.0))/1000.0, gx, gy, gz])
        time.sleep(0.0001)

    Results = np.matrix(Results)

    NextFileName = '/home/pi/Code/Prospectus/Data/Accel/60kPoints'+str(j)
    
    df = pd.DataFrame(data=Results.astype('float'))
    df.to_csv(NextFileName+'a.csv', sep=',', header=False, float_format='%.5f')
    #with open(NextFileName+'.csv','wb') as f:
    #    for line in Results:
    #        np.savetxt(f, line, fmt='%.5f')
    #Results.tofile('/home/pi/Code/Prospectus/Data/Accel/60kPoints'+str(j)+'.txt', sep=',')
    print(j, (time.time()-start) / 60.0)
    j+=1

    zip_file = zipfile.ZipFile(NextFileName+'a.zip', 'w')
    zip_file.write(NextFileName+'a.csv', compress_type=zipfile.ZIP_DEFLATED)
    zip_file.close()