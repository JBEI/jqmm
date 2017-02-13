#!/bin/bash

mkdir -p /opt/gams
wget https://d37drm4t2jghv5.cloudfront.net/distributions/24.5.6/linux/linux_x64_64_sfx.exe -O /opt/gams/linux_x64_64_sfx.exe
chmod ugo+x /opt/gams/linux_x64_64_sfx.exe
cd /opt/gams && ./linux_x64_64_sfx.exe && rm ./linux_x64_64_sfx.exe 
cp /root/gamsbatch /opt/gams/gams24.5_linux_x64_64_sfx/
chmod ugo+rx /opt/gams/gams24.5_linux_x64_64_sfx/gamsbatch
echo "export PATH=$PATH:/opt/gams/gams24.5_linux_x64_64_sfx" >> /root/.bashrc 
ln -s /opt/gams/gams24.5_linux_x64_64_sfx /gams

# note two manual steps to actually get GAMS working:
# copy gamslice.txt to /opt/gams/gams24.5_linux_x64_64_sfx/
# run ./gamsinst manually from inside /opt/gams/gams24.5_linux_x64_64_sfx/
