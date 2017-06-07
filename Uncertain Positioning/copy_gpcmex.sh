#!/bin/bash

rm -f ./Functions/Polygon/gpcmex.mexa64
cp /usr/local/MATLAB/R2016b/toolbox/map/map/private/gpcmex.mexa64 ./Functions/Polygon/
echo "gpcmex copied"
sudo chown sotiris:sotiris ./Functions/Polygon/gpcmex.mexa64
echo "Permissions fixed"

