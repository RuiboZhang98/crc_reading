# crc_reading
The multi-thread version of Gillespie algorithm for the main chain

# how to build
On Debian and Ubuntu, the code should then compile with

g++ hopscotch_MultiThread.cpp --std=c++11 -lgsl -lgslcblas -lm -o OUTPUTEXE -pthread

# how to run
Required parameters are

-o OUTPUTFILE : OUTPUTFILE is the string that will be the output CSV file name.

--cores N: N = number of threads.

--runs M: M = number of runs for each thread

Therefore MN trajectories will be simulated.

The most recent data is collected by

time nice ionice ./hopscotch -o multi_1e3x20runs_withp_notraj.csv --cores 20 --runs 1000

in which "hopscotch" is the name of OUTPUTEXE.
# output data
The output csv file contains:

1. line 1: keys of the table: tt means time(year), i means type-i crypts 
2. line 2-82: distribution of the waiting time for type-i crypts at year tt
2. line 83-163: average population for type-i crypts at year tt

# other comments
The uploaded data does not include the first line.


