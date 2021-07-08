# Configuartion is a bitch
# got mpu6050 working for 3.7 from:
#      https://pypi.org/project/py-imu-mpu6050/
#
# Initial Code based on:
#      https://pypi.org/project/mpu6050-raspberrypi/
#
# Had to enable i2c-1 with:
#      ideas from: https://circuitdigest.com/microcontroller-projects/mpu6050-gyro-sensor-interfacing-with-raspberry-pi/
#      sudo raspi-config
#      sudo apt-get install build-essential python-dev python-pip
#      sudo apt-get install python-smbus
#      sudo pip install RPi.GPIO

import mpu6050 as mp
from time import sleep
import np
sensor = mp.mpu6050(0x68)

for i in range(1250):
    accelerometer_data = sensor.get_all_data()
    print(accelerometer_data)
    sleep(0.01)


#print(sensor.get_offset())