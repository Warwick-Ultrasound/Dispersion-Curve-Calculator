# Dispersion Calculation

This is a set of scripts which allow you to calculate dispersion curves for Lamb waves in flat plates. It's been developed over many years, by many people, with each contributing to it's functionality, reliability, or speed in some way. This version was compiled by Luke in June 2023.

## Order of Operations

This section should give you an idea of what scripts to run and what you need to set at each stage. 

### cp_dispersion

Run this script first - it's the one that does the main calculation of the phase velocity vs frequency\*thickness product. You'll need to set the inputs at the top - set c_long and c_shear as te longitudinal and shear wave velocities, set the thickness of the plate, and the frequency and phase velocity step sizes and ranges. 

This script will likely take a while to run, and it'll be better on machines with a lot of cores, for example Louis.

The output of this script is a csv file with the name you gave. Inside, the first column will be frequency\*thickness, then all of the symmetric modes, then the antisymmetric modes. On completion of the script, "SPLIT = N" will be printed in the console - note this dfown as it's the column number that separates the symmetric and antisymmetric modes, and you'll need it in a bit.

### jump_remover

Almost always, the output of cp_dispersion will have jumps in the data, and msising postions of some lines where NaNs are in the array. You need to fix this before you move onto the next step, so jump_remover does the fixing. It'll ouput another file containing all the modes with the same name but "GAP_FILLED" appended to the end.

Check that it looks sensible before you continue.

### calc_Cg

Once you've filled the gaps just run calc_Cg with the GAP_FILLED filename inputted at the top of the script. This should just work but you might have to play around with the filtering a bit if there are lots of wiggly bits in the data.

If you're struggling to get to high enough frequency\*thickness, one thing I found that helps is to calculate quite a bit past the frequency thickness you need right from the start, then crop of the unstable bit at the end.