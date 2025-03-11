#!/bin/zsh
 
start=-1.5
end=1.5
step=0.05

OUTFILE="energy_XXZ.txt"

echo -e "Delta\t\tEnergy" > "$OUTFILE"

for Delta in $(seq $start $step $end); do
    sed -i '' "s/^ *Delta *=.*/Delta=$Delta/" parameter.py
    echo "$Delta"
    /Users/swarnamoyghosh/anaconda3/envs/drudge/bin/python3 XXZ_energy.py
    
done

echo "completed!"
