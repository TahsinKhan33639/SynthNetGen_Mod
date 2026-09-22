#!/usr/bin/env bash
set -e

cd -- "$(dirname -- "$(realpath -- "${BASH_SOURCE[0]}")")"

g++ -O2 cleanup.cpp -o cleanup
g++ -O2 expfy.cpp -o expfy
g++ -O2 onb.cpp -o onb
g++ -O2 bnb.cpp -o bnb
g++ -O2 gen_deg2.cpp -o gen_deg2

g++ -O2 bxk4f_rseas.cpp -o bxk4f_rseas
g++ -O2 bxk4f_rmse.cpp -o bxk4f_rmse
g++ -O2 bxk4one_rseas.cpp -o bxk4one_rseas
g++ -O2 bxk4one_rmse.cpp -o bxk4one_rmse

g++ -O2 otb.cpp -o otb2
g++ -O2 orca.cpp -o orca
g++ -O2 evs.cpp -o evs
g++ -O2 evc1.cpp -o evc1
g++ -O2 evc2.cpp -o evc2
g++ -O2 evc3.cpp -o evc3

chmod +x rpll.sh rpllb.sh
