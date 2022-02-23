# leap_sims

### How to run leap_sims_detector_test
1.  set the environment:  
  `source /cvmfs/sft.cern.ch/lcg/views/LCG_98python3/x86_64-centos7-gcc10-opt/setup.sh`  

2.  make a build directory and change to the build directory:  
  `mkdir build`  
  `cd build`  

3. build the application:  
  `cmake ../`  
  `make`  
4. run leap_sims with a .mac file
  `./leap_sims  -m macfile.mac -f outFileName -t outType -v version -b <on/off>`
  for example use: `./leap_sims  -m test45.mac -f test_result_file.root -t bunch -v Cal -b on`

  `-m` specifies the name of the macrofile and starts the run in batch mode. At the moment there a two files you can use `test45.mac` and `test_init_vis.mac` that          is a copy of `init_vis.mac` without starting the visualization.
  `-f` specifies the name of the ttree saved. At that stage you have to give a different name for every run because a .gdml file will be created to check the              geometry. 
  `-t` specifies the what is saved to the ttree. Use `bunch` at the moment! (bunch is default)   
  `-v` specifies the version of the code that is used. `PolCal` simulates polarimeter and calorimeter. `Pol` simulates only the polarimeter. `Cal` simulates only          the calorimeter.
  `-b` you can switch the beam line `on` and `off`

  **Important**: At the moment just version type `Cal` can be used!

5. run visualization
  to start with visualization just don't supply the .mac file
  `./leap_sims -f outFileName -t outType -v version -b <on/off>`
  for examlple use: `./leap_sims -f test_result_file.root -t bunch -v Cal -b on`
