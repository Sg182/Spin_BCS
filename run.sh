#!/bin/bash
 
start=-2.0
end=2.0
step=0.1

for Delta in $(seq $start $step $end); do
    /Users/swarnamoyghosh/anaconda3/bin/python XXZ_energy.py $Delta
done

echo "Done!"
