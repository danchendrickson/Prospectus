sudo apt-get update
sudo apt full-upgrade -y
sudo apt-get install build-essential python-dev python-pip -y
sudo apt-get install code -y
sudo apt-get install python-smbus -y

sudo pip install smbus2 #i2c-tools
#sudo pip install py-imu-mpu6050
sudo pip install numpy
sudo pip install matplotlib
sudo pip install RPi.GPIO 

sudo raspi-config
