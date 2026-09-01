#!/usr/bin/env bash
set -e

g++ cleanup.cpp -o cleanup
g++ expfy.cpp -o expfy
g++ onb.cpp -o bnb
g++ gen_deg2.cpp -o gen_deg2

g++ bxk4f_rseas.cpp -o bxk4f_rseas
g++ bxk4f_rmse.cpp -o bxk4f_rmse
g++ bxk4one_rseas.cpp -o bxk4one_rseas
g++ bxk4one_rmse.cpp -o bxk4one_rmse

g++ otb.cpp -o otb2
g++ orca.cpp -o orca
g++ evs.cpp -o evs
g++ evc1.cpp -o evc1
g++ evc2.cpp -o evc2

chmod +x rpll.sh
