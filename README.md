To make the SGP4 model easier to use, this project added the topology generation module, which can get the information of User Terminals and satellites by reading TXT files.


step 1:
In example/passpredict.cc, specify the locations of user_terminals.txt and tles.txt
Step 2:
mkdir bulid, cmake .. make
step 3:
./predict