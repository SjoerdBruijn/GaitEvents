
Code to calculate heelstrike and toe offs, and correct these using a simple GUI. The first option is to calculate the heelstrikes (and toe-off) and check them based on ```cop``` data as in the algorithm from Roerdink[1], as explained in this video (the GUI is also explained there):

<iframe width="500" height="300" src="https://www.youtube.com/embed/LWktHNFrAtM?si=0d3FJcgGZDkwfXs_" frameborder="0" allowfullscreen></iframe>

If you download the entire repo, and have your path set to the repo, you should be able to run the following code examples; in Python:
(for some reason, in this example, you will need to recalculate the evenst with the win set to 0.1 in the GUI, which is not the case in the Matlab example; working to fix this)

Note; for the CoP, the columns should be: 1.ML(+ = right) the  2.AP(+ = forward)
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

in Matlab:
```
clear all; close all;clc
addpath(genpath('GaitEventsMatlab'))
filename= 'ExampleData/ExampleCoP.csv'
fs      = 200;
cop     = csvread(filename);
events  = calc_events(cop, fs);
events  = check_events(events,cop,fs);
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
filename= 'ExampleData/Exampletraj.mat';
data    = load(filename);
fs      = data.fs;
traj    = data.traj;
events  = calc_events(traj, fs);
events  = check_events(events,traj,fs)
```

if you dont have the ```traj``` structure, but you do have a C3D file, here is a trick to get those data in a ```traj``` structure.
in Python (this example requires the ezc3d package, see here https://github.com/pyomeca/ezc3d). Notes; this particular C3D file contains a period of standing still first, and you will have to remove quite a few of wrongle detected events in the beginning of the file; all points before 1400 samples are basically wrong. This example also leads to some weird results; recalculate with win=2 when checking gives you good results:
```
import numpy as np
from ezc3d import c3d
import GaitEventsPython.GaitEvents as GE


filename = 'ExampleData/ExampleC3D.c3d'
c = c3d(filename)
point_data = c['data']['points']
videoframeRate =c['header']['points']['frame_rate']
LHI = c['parameters']['POINT']['LABELS']['value'].index('LHEE') #actually shoudl contain try statement, but I know this label is there, so lazy...
RHI = c['parameters']['POINT']['LABELS']['value'].index('RHEE')

dt = np.dtype([('name', 'U20'), ('joint_name', 'U20'), ('origin', 'O')])
traj=dict()
traj['segment']=np.zeros((1,2), dtype=dt)
traj['segment'][0,0]['name'] = 'Left foot';
traj['segment'][0,1]['name'] = 'Right foot';
traj['segment'][0,0]['joint_name'] = '[]';
traj['segment'][0,1]['joint_name'] = '[]';
traj['segment'][0,0]['origin']  = np.expand_dims(point_data[0:3,LHI,:],axis=1)
traj['segment'][0,1]['origin']  = np.expand_dims(point_data[0:3,RHI,:],axis=1)
events,_= GE.calc_events(traj, videoframeRate)
checker = GE.GaitEventChecker(events,traj,videoframeRate)
events  = checker.events
```` 
 
in Matlab (you will see same problems due to stanidng still, however, not the win problem): 

```
clear all; close all;clc
addpath(genpath('GaitEventsMatlab'))
filename = 'ExampleData/ExampleC3D.c3d'
[POINTdat,VideoFrameRate,ANALOGdat,AnalogFrameRate,Event,ParameterGroup,CameraInfo,ResidualError]...
= readC3D_mhs(filename);
traj.segment(1).name{1} = 'Left foot';
traj.segment(2).name{1} = 'Right foot';
traj.segment(1).joint_name = '';
traj.segment(2).joint_name = '';
traj.segment(1).origin  = GetMarkerDataVicon(POINTdat,ParameterGroup,{'LHEE'});
traj.segment(2).origin  = GetMarkerDataVicon(POINTdat,ParameterGroup,{'RHEE'});
events  = calc_events(traj, VideoFrameRate);
events  = check_events(events,traj,VideoFrameRate);

```



Lastly, you can also use this in several other ways, which just use the GUI in a smart way. For instance, if the above method of calculating heelstrikes based on the vertical position of the foot markers doesnt work, you can also calculate events based on the AP position of the foot markers (or some other thinsg). Again, Matlab only for now (but should be easy to implement in Python as well)
Note, same file as above, so same problem... here we solve that by simply ignoring the first part of the timeseries
Also NOTE: you CAN NOT YET use the "recalculate" button if you use a trick like this. 
in Python:
```
import numpy as np
from scipy.signal import find_peaks
from ezc3d import c3d
import GaitEventsPython.GaitEvents as GE


filename = 'ExampleData/ExampleC3D.c3d'
c = c3d(filename)
point_data = c['data']['points']
videoframeRate =c['header']['points']['frame_rate']
LHI = c['parameters']['POINT']['LABELS']['value'].index('LHEE') #actually shoudl contain try statement, but I know this label is there, so lazy...
RHI = c['parameters']['POINT']['LABELS']['value'].index('RHEE')
Lfoot = point_data[1,LHI,2999:].T
Rfoot = point_data[1,RHI,2999:].T
events=dict()
events['lhs'],_  = find_peaks(-Lfoot,distance=videoframeRate)# we use max forward position of the foot as heelstrike
events['lto'],_  = find_peaks(Lfoot,distance=videoframeRate)# and max backwoard position of the foot as toe off
events['rhs'],_  = find_peaks(-Rfoot,distance=videoframeRate) # we use max forward position of the foot as heelstrike
events['rto'],_ = find_peaks(Rfoot,distance=videoframeRate) # and max backwoard position of the foot as toe off
checker = GE.GaitEventChecker(events,np.vstack((Lfoot, Rfoot)).T,videoframeRate)
events  = checker.events
````
In Matlab
```
clear all; close all; clc
addpath(genpath('GaitEventsMatlab'))
filename = 'ExampleData/ExampleC3D.c3d'

[POINTdat,VideoFrameRate,ANALOGdat,AnalogFrameRate,Event,ParameterGroup,CameraInfo,ResidualError]...
    = readC3D_mhs(filename);

Lfoot   = GetMarkerDataVicon(POINTdat,ParameterGroup,{'LHEE'});
Rfoot   = GetMarkerDataVicon(POINTdat,ParameterGroup,{'RHEE'});
Lfoot   = -squeeze(Lfoot(2,:,3000:end)); % we only use x position, and x was negatively forward here. 
Rfoot   = -squeeze( Rfoot(2,:,3000:end));
[Lmaxout, Lminout] = peakfind(Lfoot,VideoFrameRate);
events.lhs  = Lmaxout(:,1); % we use max forward position of the foot as heelstrike
events.lto  = Lminout(:,1); % and max backwoard position of the foot as toe off
[Rmaxout, Rminout] = peakfind(Rfoot,VideoFrameRate);
events.rhs  = Rmaxout(:,1); % we use max forward position of the foot as heelstrike
events.rto  = Rminout(:,1); % and max backwoard position of the foot as toe off
events  = check_events(events,[Lfoot,Rfoot],VideoFrameRate);



```
Or you could do it with ground reaction forces, and some kind of % of GRF; 
In python:
```
import GaitEventsPython.GaitEvents as GE
import numpy as np

# Load data from CSV
filename = 'ExampleData/ExampleForces.csv'
fs = 1000
threshold_force = 100  # or some other number, just what you like

# Read CSV data into a NumPy array using genfromtxt
forces = np.genfromtxt(filename, delimiter=',')

# Extract relevant columns (indexing starts from 0 in Python)
Fzleft = forces[:, 2]
Fzright = forces[:, 5]

# Calculate events using numpy
events = {
    'lhs': np.where(np.diff((Fzleft > threshold_force).astype(int)) == 1)[0],  # left heelstrike
    'lto': np.where(np.diff((Fzleft > threshold_force).astype(int))  == -1)[0], # left toe-off
    'rhs': np.where(np.diff((Fzright > threshold_force).astype(int))  == 1)[0], # right heelstrike
    'rto': np.where(np.diff((Fzright > threshold_force).astype(int))  == -1)[0] # right toe-off
}
checker = GE.GaitEventChecker(events,forces[:,(2 ,5)],fs)
events  = checker.events


```
in Matlab: 
```
clear all; close all;clc
addpath(genpath('GaitEventsMatlab'))
filename= 'ExampleData/ExampleForces.csv'
fs      = 1000
treshold_force = 100; % or some other number, just what you like
forces  = csvread(filename);
Fzleft  = forces(:,3);
Fzright = forces(:,6);
events.lhs = find(diff(Fzleft>treshold_force)==1); % left heelstrike is when the force goes over treshold
events.lto = find(diff(Fzleft>treshold_force)==-1);
events.rhs = find(diff(Fzright>treshold_force)==1); % left heelstrike is when the force goes over treshold
events.rto = find(diff(Fzright>treshold_force)==-1);
events  = check_events(events,[Fzleft,Fzright],fs)
```

references:
[1] Roerdink
[2] pijnappels
