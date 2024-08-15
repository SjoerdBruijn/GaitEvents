
Code to calculate heelstrike and toe offs, and correct these using a simple GUI. The first option is to calculate the heelstrikes (and toe-off) and check them based on ```cop``` data as in the algorithm from Roerdink[1], as explained in this video (the GUI is also explained there):

<iframe width="500" height="300" src="https://www.youtube.com/embed/LWktHNFrAtM?si=0d3FJcgGZDkwfXs_" frameborder="0" allowfullscreen></iframe>

If you download the entire repo, and have your path set to the repo, you should be able to run the following code examples; in Python:
```
import GaitEventsPython.GaitEvents as GE
import numpy as np
filename= 'ExampleData/ExampleCoP.csv'
cop     = np.genfromtxt(filename, delimiter=",")
fs      = 200
events,_  = GE.calc_events(cop, fs)
checker = GE.GaitEventChecker(events,cop,fs)
events  = checker.events
```
(for some reason, in the above example, you will need to recalculate the evenst with the win set to 0.1, which is not the case in the Matlab example; working to fix this)

in Matlab:
```
clear all; close all;clc
addpath(genpath('GaitEventsMatlab'))
filename= 'ExampleData/ExampleCoP.csv'
fs      = 200
cop     = csvread(filename)
events  = calc_events(cop, fs)
events  = check_events(events,cop,fs)
```

You could also calculate heelstrikes (and toe offs) and check them based on a VU 3d model ```traj``` data structure, in which case the vertical position and velocity algorithm will be used; check the video here; 
<iframe width="500" height="300" src="https://www.youtube.com/embed/WD3iCJZIdJo?si=b3kCsuAG8ZeO0sfS" frameborder="0" allowfullscreen></iframe>

in Python:
```
import GaitEventsPython.GaitEvents as GE
from scipy.io import loadmat
filename= 'ExampleData/Exampletraj.mat'
data    = loadmat(filename)
traj    = data["traj"][0][0]
fs      = int(data["fs"][0][0])
events,_= GE.calc_events(traj, fs)
checker = GE.GaitEventChecker(events,traj,fs)
events  = checker.events
```

in Matlab:
```
clear all; close all;clc
addpath(genpath('GaitEventsMatlab'))
filename= 'ExampleData/Exampletraj.mat'
data    = load(filename)
fs      = data.fs;
traj    = data.traj;
events  = calc_events(traj, fs);
events  = check_events(events,traj,fs)
```

if you dont have the ```traj``` structure, but you do have a C3D file, here is a trick to get those data in a ```traj``` structure in Matlab (sorry, no time to do this in Python yet, but should not be too hard to implement) 

```
clear all; close all;clc


filename = 'ExampleData/ExampleC3D.c3d'
[POINTdat,VideoFrameRate,ANALOGdat,AnalogFrameRate,Event,ParameterGroup,CameraInfo,ResidualError]...
= readC3D_mhs(filename);
traj.segment(1).name{1} = 'Left foot';
traj.segment(2).name{1} = 'Right foot';
traj.segment(1).joint_name = '';
traj.segment(2).joint_name = '';
traj.segment(1).origin  = GetMarkerDataVicon(POINTdat,ParameterGroup,{'LHEE'});
traj.segment(2).origin  = GetMarkerDataVicon(POINTdat,ParameterGroup,{'RHEE'});
events  = calc_events(traj, VideoFrameRate)
events  = check_events(events,traj,VideoFrameRate)

```
(NOTE: the above file contains a period of standing still first, and you will have to remove quite a few of wrongle detected events in the beginning of the file; all points before 1400 samples are basically wrong. )


Lastly, you can also use this in several other ways, which just use the GUI in a smart way. For instance, if the above method of calculating heelstrikes based on the vertical position of the foot markers doesnt work, you can also calculate events based on the AP position of the foot markers (or some other thinsg). Again, Matlab only for now (but should be easy to implement in Python as well)
NOTE: same file as above, so similar problem. Also NOTE: you CAN NOT YET use the "recalculate" button if you use a trick like this. 
```
clear all; close all; clc
filename = 'ExampleData/ExampleC3D.c3d'

[POINTdat,VideoFrameRate,ANALOGdat,AnalogFrameRate,Event,ParameterGroup,CameraInfo,ResidualError]...
    = readC3D_mhs(filename);

Lfoot   = GetMarkerDataVicon(POINTdat,ParameterGroup,{'LHEE'});
Rfoot   = GetMarkerDataVicon(POINTdat,ParameterGroup,{'RHEE'});
Lfoot   = squeeze(Lfoot(2,:,:)); % we only use x position.
Rfoot   = squeeze( Rfoot(2,:,:));
[Lmaxout, Lminout] = peakfind(Lfoot,VideoFrameRate);
events.lhs  = Lmaxout(:,1); % we use max forward position of the foot as heelstrike
events.lto  = Lminout(:,1); % and max backwoard position of the foot as toe off
[Rmaxout, Rminout] = peakfind(Rfoot,VideoFrameRate);
events.rhs  = Rmaxout(:,1); % we use max forward position of the foot as heelstrike
events.rto  = Rminout(:,1); % and max backwoard position of the foot as toe off
events  = check_events(events,[Lfoot,Rfoot],VideoFrameRate)



```
Or you could do it with ground reaction forces, and some kind of % of GRF; 
Note; this is not the best example dataset, as it contains quite some steps which are on both belts.   

```
clear all; close all;clc
addpath(genpath('GaitEventsMatlab'))
filename= 'ExampleData/ExampleForces.csv'
fs      = 200
treshold_force = 100; % or some other number, just what you like
forces  = csvread(filename);
Fzleft  = forces(:,3);
Fzright = forces(:,6);
events.lhs = find(diff(Fzleft>treshold_force)==1); % left heelstrike is when the force goes over treshold
events.lto = find(diff(Fzleft>treshold_force)==-1);
events.rhs = find(diff(Fzright>treshold_force)==1); % left heelstrike is when the force goes over treshold
events.rto = find(diff(Fzright>treshold_force)==-1);
events  = check_events(events,[Fzleft,Fzright],fs)


references:
[1] Roerdink
[2] pijnappels