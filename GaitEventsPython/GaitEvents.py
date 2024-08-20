import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy import signal as sig
from scipy.signal import find_peaks
from scipy import interpolate



class GaitEventChecker:
    def __init__(self, events, signal, fs, maxsubplots=2):

        self.rawsignal = signal
        self.fs = fs
        _,self.signal = calc_events(self.rawsignal, self.fs)
        try:
            self.events,_ = order_events(events)
        except:
            self.events = events
        self.maxsubplots = maxsubplots
        self.notgood = True
        self.is_done = False
        self.i_plot = 0
        self.ws = 10  # window size in seconds
        self.maxwindow = self.ws * self.fs

        if self.maxwindow > len(self.signal):
            self.ws = round(len(self.signal) / self.fs / 2)
            self.maxwindow = self.ws * self.fs

        self.maxplots = int(np.ceil(len(self.signal) / self.maxwindow)) + 1

        # Initialize the Tkinter window
        self.root = tk.Tk()
        self.root.title("Gait Event Checker")

        self.create_widgets()
        self.update_plot()
        self.root.protocol("WM_DELETE_WINDOW", self.exit_no_good)
        self.root.mainloop()

    def create_widgets(self):
        # Create the main frame
        self.main_frame = ttk.Frame(self.root)
        self.main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Create a frame for the plots and the side buttons
        plot_and_buttons_frame = ttk.Frame(self.main_frame)
        plot_and_buttons_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Create the side frame for the first set of buttons
        side_button_frame = ttk.Frame(plot_and_buttons_frame)
        side_button_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)

        # First set of buttons next to the figure
        side_buttons = [
            ("Done, exit", self.done),
            ("No good data", self.exit_no_good),
            ("Calc Events", self.recalculate_events),
            ("OK, next window", self.next_window),
            ("One window back", self.prev_window),
            ("Delete selected point(d)", self.delete_point),
            ("Add selected point(a)", self.add_point)]

        for text, command in side_buttons:
            button = ttk.Button(side_button_frame, text=text, command=command)
            button.pack(fill=tk.X, pady=5)

        # Adding fc_label and win_label entries
        self.fc_label = ttk.Label(side_button_frame, text="fc:")
        self.fc_label.pack(anchor='w', pady=2)

        self.fc_entry = ttk.Entry(side_button_frame)
        self.fc_entry.pack(fill=tk.X, pady=2)
        # self.fc_entry.insert(0,None)
       

        self.win_label = ttk.Label(side_button_frame, text="win:")
        self.win_label.pack(anchor='w', pady=2)

        self.win_entry = ttk.Entry(side_button_frame)
        self.win_entry.pack(fill=tk.X, pady=2)
        self.win_entry.insert(0,1)

        # Create a frame for the figure
        figure_frame = ttk.Frame(plot_and_buttons_frame)
        figure_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Initialize matplotlib figure and two subplots
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax1 = self.fig.add_subplot(211)  # Top subplot
        self.ax2 = self.fig.add_subplot(212)  # Bottom subplot
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('key_release_event', self.KeyDownListener)

        self.canvas = FigureCanvasTkAgg(self.fig, master=figure_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Create the bottom frame for the second set of buttons
        bottom_frame = ttk.Frame(self.main_frame)
        bottom_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=10)

        # Second set of buttons arranged in two rows and four columns
        bottom_buttons = [
            ("Insert lhs", lambda: self.dummy_function("Insert lhs"),"magenta"),
            ("Insert rto", lambda: self.dummy_function("Insert rto"),"green"),
            ("Insert rhs", lambda: self.dummy_function("Insert rhs"),"black"),
            ("Insert lto", lambda: self.dummy_function("Insert lto"),"blue"),
            ("Delete lhs", lambda: self.dummy_function("Delete lhs"),"magenta"),
            ("Delete rto", lambda: self.dummy_function("Delete rto"),"green"),
            ("Delete rhs", lambda: self.dummy_function("Delete rhs"),"black"),
            ("Delete lto", lambda: self.dummy_function("Delete lto"),"blue"),
        ]

        for i, (text, command,color) in enumerate(bottom_buttons):
            button = tk.Button(bottom_frame, text=text, command=command, bg=color,fg=color)
            button.grid(row=i // 4, column=i % 4, padx=5, pady=5, sticky="ew")

        # Configure the grid to expand
        for i in range(4):
            bottom_frame.grid_columnconfigure(i, weight=1)
            
    def on_click(self, event):
        if event.inaxes == self.ax1:
            self.i_plot = round(
                self.maxplots * (event.xdata / self.ax1.get_xlim()[1]))
            self.update_plot()
    
    def KeyDownListener(self,event):
        if event.key == 'd':
            self.delete_point()
        elif event.key == 'a':
            self.add_point()
        elif event.key == 'right':
            self.next_window()
        elif event.key == 'left':
            self.prev_window()
        elif event.key == 'up':
            self.ratP = self.i_plot / self.maxplots
            if self.ws < int(len(self.signal) / 2) - 1:
                self.ws += 1
                self.maxwindow = self.ws * self.fs
                self.maxplots = int(np.ceil(len(self.signal) / self.maxwindow) + 1)
            self.i_plot = int(self.ratP * self.maxplots)
            self.update_plot()
            
        elif event.key == 'down':
            self.ratP = self.i_plot /self. maxplots
            if self.ws > 1:
                self.ws -= 1
                self.maxwindow = self.ws * self.fs
                self.maxplots = int(np.ceil(len(self.signal) / self.maxwindow) + 1)
            self.i_plot = int(self.ratP * self.maxplots)
            self.update_plot()
            
    def update_plot(self):
         self.ax1.clear()
         self.ax1.plot(np.linspace(1, max(np.shape(self.signal)), len(
             np.diff(self.events['lto']))), np.diff(self.events['lto']), color='blue')
         self.ax1.plot(np.linspace(1, max(np.shape(self.signal)), len(
             np.diff(self.events['lhs']))), np.diff(self.events['lhs']), color='magenta')
         self.ax1.plot(np.linspace(1, max(np.shape(self.signal)), len(
             np.diff(self.events['rto']))), np.diff(self.events['rto']), color='green')
         self.ax1.plot(np.linspace(1, max(np.shape(self.signal)), len(
             np.diff(self.events['rhs']))), np.diff(self.events['rhs']), color='black')
        
         self.ax1.axhline(y=np.nanmean(np.diff(self.events['lhs'])) -
                             np.std(np.diff(self.events['lhs'])), color='gray', linestyle='--')
         self.ax1.axhline(y=np.nanmean(np.diff(self.events['lhs'])) +
                             np.std(np.diff(self.events['lhs'])), color='gray', linestyle='--')
        
         Meve = np.column_stack((np.concatenate((self.events['lhs'], self.events['rto'], self.events['rhs'], self.events['lto'])),
                                 np.concatenate((np.zeros_like(self.events['lhs']), np.ones_like(self.events['rto']),
                                                 np.ones_like(self.events['rhs']) * 2, np.ones_like(self.events['lto']) * 3))))
         idxwrong = []
         if len(Meve) > 32:
             Meve = Meve[Meve[:, 0].argsort()]
             Morder = np.column_stack((np.arange(1, len(Meve) + 1) % 4,
                                       np.arange(2, len(Meve) + 2) % 4,
                                       np.arange(3, len(Meve) + 3) % 4,
                                       np.arange(4, len(Meve) + 4) % 4))
             dMeve = Morder - Meve[:, 1][:, np.newaxis]
             eveordercol = np.argmin(np.sum(dMeve != 0, axis=0))
             ieve_wrong = np.argmax(
                 np.abs(Morder[:, eveordercol] - Meve[:, 1]) > 0)
             idxwrong = Meve[ieve_wrong, 0]
             evewrong = Meve[ieve_wrong, 1] + 1
             self.ax1.axvline(x=idxwrong, color='red',
                                 linestyle=':', linewidth=2)
        
         lijnkleur = ['black', 'blue', 'black', 'blue', 'cyan', 'black']
        
         domein = np.arange(round((self.i_plot * self.maxwindow)),
                            round((self.i_plot + 1) * self.maxwindow))  
        
         if domein[0] < 0:
             domein = domein + abs(domein[0]) + 1
         if domein[-1] > max(np.shape(self.signal)):
             domein = domein - (domein[-1] - max(np.shape(self.signal)))-1
             
         self.ax1.set_ylim(np.nanmean(np.diff(
             self.events['lhs'])) - 0.5*self.fs, np.nanmean(np.diff(self.events['lhs'])) + 0.5*self.fs)
         self.ax1.axvline(x=domein[0], color='black', linewidth=2)
         self.ax1.axvline(x=domein[-1], color='black', linewidth=2)
         self.ax2.clear()
         for j in range(1, self.maxsubplots):
             domein_j = domein[
                 round(len(domein) * (1 / (self.maxsubplots - 1)) * (j - 1)):round(len(domein) * (j / (self.maxsubplots - 1)))]
        
             i_lhs = np.where((self.events['lhs'] > domein_j[0]) & (
                 self.events['lhs'] < domein_j[-1]))[0]
             i_rto = np.where((self.events['rto'] > domein_j[0]) & (
                 self.events['rto'] < domein_j[-1]))[0]
             i_rhs = np.where((self.events['rhs'] > domein_j[0]) & (
                 self.events['rhs'] < domein_j[-1]))[0]
             i_lto = np.where((self.events['lto'] > domein_j[0]) & (
                 self.events['lto'] < domein_j[-1]))[0]
        
             # plt.subplot(self.maxsubplots, 1, j + 1)
             for k in range(self.signal.shape[1]):
                 self.ax2.plot(
                     domein_j, self.signal[domein_j, k], color=lijnkleur[k])
                 self.ax2.plot(self.events['lto'][i_lto],
                                  self.signal[self.events['lto'][i_lto], k], 'o', color='blue')
                 self.ax2.plot(self.events['lhs'][i_lhs],
                                  self.signal[self.events['lhs'][i_lhs], k], 'o', color='magenta')
                 self.ax2.plot(self.events['rto'][i_rto],
                                  self.signal[self.events['rto'][i_rto], k], 'o', color='green')
                 self.ax2.plot(self.events['rhs'][i_rhs],
                                  self.signal[self.events['rhs'][i_rhs], k], 'o', color='black')
            
             if idxwrong and domein_j[0] < idxwrong < domein_j[-1]:
                 self.ax2.plot(idxwrong, self.signal[idxwrong, 0],
                                  'o', color='red', markersize=10)
                 self.ax2.plot(idxwrong, self.signal[idxwrong, 1],
                                  'o', color='red', markersize=10)
        
             self.ax2.set_title('points check correct using options left')
             self.ax2.set_xlim(domein_j[0], domein_j[-1])
             
         if not (len(self.events['lto']) == 
                 len(self.events['rto']) == 
                 len(self.events['rhs']) == 
                 len(self.events['lhs'])):
             self.ax1.set_title(
                 "Not all points are in the order lhs=>rto=>rhs=>lto, keep correcting") 
         elif (any(self.events['rto'] - self.events['lhs'] < 0) |
               any(self.events['rhs'] - self.events['rto'] < 0) |
               any(self.events['lto'] - self.events['rhs'] < 0) |
               any(self.events['lhs'][1:] - self.events['lto'][:-1] < 0)):
             self.ax1.set_title(
                 "Not all points are in the order lhs=>rto=>rhs=>lto, keep correcting")
         else:
             self.ax1.set_title(
                 "All points are in the order lhs=>rto=>rhs=>lto :)")

         self.canvas.draw()

    def done(self):
        self.is_done = True
        plt.close(self.fig)
        self.root.quit()
        self.root.destroy()

    def exit_no_good(self):
        self.notgood = True
        for key in self.events:
            self.events[key] = []
        plt.close(self.fig)
        self.root.quit()
        self.root.destroy()

    def recalculate_events(self):
        # Placeholder for recalculating events
        print("Recalculating events...")
        if self.fc_entry.get() =='':
            fc = None
        else:
            fc = float(self.fc_entry.get())
        
        tmp = float(self.win_entry.get()) 
        self.events,self.signal=calc_events(self.rawsignal, self.fs, tmp,fc)
        self.update_plot()

    def next_window(self):
        if self.i_plot < self.maxplots - 1:
            self.i_plot += 1
            self.update_plot()

    def prev_window(self):
        if self.i_plot > 0:
            self.i_plot -= 1
            self.update_plot()

    def delete_point(self):
        print("Deleting selected point...")
        self.ax1.set_title("Click on the plot to select a point to delete.")
        tmp = self.fig.ginput(1, show_clicks=True, timeout=-1)
        x=tmp[0][0]
        distances = {
             'lto': np.abs(self.events['lto'] - x),
             'rto': np.abs(self.events['rto'] - x),
             'lhs': np.abs(self.events['lhs'] - x),
             'rhs': np.abs(self.events['rhs'] - x),
         }
        
        # Find the minimum distance and corresponding event category
        min_dist = {key: np.min(distances[key]) if len(
             distances[key]) > 0 else np.inf for key in distances}
        min_idx = {key: np.argmin(distances[key]) if len(
             distances[key]) > 0 else -1 for key in distances}
        variable = min(min_dist, key=min_dist.get)
        
        if min_idx[variable] != -1:
            self.events[variable] = np.delete(
                self.events[variable], min_idx[variable])
        self.update_plot()
         
    def add_point(self):
        # Placeholder for adding a selected point
        print("Adding selected point...")
        self.ax1.set_title("Click on the plot to add a point.")
        tmp = self.fig.ginput(1, show_clicks=True, timeout=-1)

        x=tmp[0][0]
        y=tmp[0][1]
        def mindist(target, data):
            """
            Find the minimal distance to the target.
            
            Parameters:
                target (float or np.ndarray): Target value from where distance is computed (scalar or vector).
                data (np.ndarray): The data from which the distance is measured (matrix or vector).
                
            Returns:
                max_value (float or np.ndarray): Value of data point closest to the target.
                index (int or np.ndarray): Index of closest point in data (per column if matrix).
            """
            # Replicate the target to match the size of the data
            target = np.tile(target, data.shape)
            
            # Compute the minimum distance and the corresponding index
            abs_diff = np.abs(target - data)
            index = np.argmin(abs_diff, axis=0)
            
            if len(data.shape) > 1 and data.shape[0] > 1 and data.shape[1] > 1:
                # If the data is a 2D array (matrix), return the diagonal of the selected points
                max_value = data[index, np.arange(data.shape[1])]
            else:
                # If the data is 1D, or a degenerate case, return the selected point
                max_value = data[index]
            
            return max_value, index

        # Find closest event-point to input x
        _, e1_1 = mindist(x, self.events['lhs'])
        _, e1_2 = mindist(x, self.events['rto'])
        _, e1_3 = mindist(x, self.events['rhs'])
        _, e1_4 = mindist(x, self.events['lto'])
        
        event_values = np.array([self.events['lhs'][e1_1], self.events['rto'][e1_2], self.events['rhs'][e1_3], self.events['lto'][e1_4]])
        val_e, e2 = mindist(x, event_values)
        
        # Determine which variable to use
        if x - val_e > 0:
            vars_list = ['rto', 'rhs', 'lto', 'lhs']
        else:
            vars_list = ['lto', 'lhs', 'rto', 'rhs']
        
        # Point to add variable
        variable = vars_list[e2]
        
        # Round x to calculate minimal distance to top
        x = round(x)
        
        # Find signal closest to selected point
        _, sign = mindist(y, self.signal[x, :])
        
        # Let the peak find window size be a function of part of the signal shown
        R = round(self.maxwindow / 20)
        
        # Don't survey complete signal, only 1 second surrounding selected point
        range_vals = np.arange(max(0, x - R), min(len(self.signal), x + R))
        
        # Select the window for the peak find
        Peak_ws = round(0.1 * len(range_vals))
        
        # Find closest peak to the selected point
        peaks_max, _ = find_peaks(self.signal[range_vals, sign], distance=Peak_ws)
        peaks_min, _ = find_peaks(-self.signal[range_vals, sign], distance=Peak_ws)
        
        if len(peaks_max) > 0 and len(peaks_min) > 0:
            all_peaks = np.concatenate((range_vals[peaks_max], range_vals[peaks_min]))
            _, closest_peak_idx = mindist(x, all_peaks)
            new_point = all_peaks[closest_peak_idx]
        elif len(peaks_max) == 0 and len(peaks_min) > 0:
            new_point = range_vals[peaks_min.min()]
        elif len(peaks_min) == 0 and len(peaks_max) > 0:
            new_point = range_vals[peaks_max.max()]
        else:
            # If none is found, use the selected point as new point
            new_point = x
        
        # Add new point to the existing points
        self.events[variable] = np.sort(np.append(self.events[variable], new_point))
        self.update_plot()

    def dummy_function(self, button_name):
        # This function now prints which button was pressed
        print(f"{button_name} now")
        variable=button_name[-3:].lower()
    
        # self.axs[0].set_title("Click on the plot to add a point.")
        tmp = self.fig.ginput(1, show_clicks=True, timeout=-1)    
        new_point = tmp[0][0]
        if button_name[0:3]=="Ins":
            self.events[variable] = np.sort(np.append(self.events[variable], int(new_point)))
        else:
            ind=np.where(new_point<self.events[variable])
            self.events[variable] = np.delete(
                self.events[variable], ind[0][0])
        self.update_plot()
        
def calc_events(signal, fs, tmpcor=1,fc=None):
    """
    Function to calculate heel strikes and toe offs

    Parameters:
    signal (numpy.ndarray or dict): Center of Pressure data in columns (or traj structure)
    fs (float): Sample frequency of your data
    tmpcor (float, optional): Correction for the estimation of the step length. Defaults to 1.

    Returns:
    dict: Structure containing the indices of the gait events of the CoP data
    """
    events = {}

      
    if len(signal)==1 or len(signal)==2:
        traj = signal
        signal = np.empty((len(traj["segment"][0,0]['origin'][2, 0, :]),))
        for i_seg in range(0,len(traj["segment"][0,])):
            if 'foot' in str(traj["segment"][0,i_seg]["name"]).lower():
                if not str(traj["segment"][0,i_seg]["joint_name"])=='[]':
                    for i_joint in range(0,len(traj["segment"][0,i_seg]["joint_name"][0][:])):
                        flag=0
                        if 'heel' in str(traj["segment"][0,i_seg]["joint_name"][0][i_joint]).lower():
                            z = traj["segment"][0,i_seg]['joint'][2,i_joint,:]
                            flag=1
                    if flag==0:
                        z = traj["segment"][0,i_seg]['joint'][2,2,:]

                        
                else:
                    z = traj["segment"][0,i_seg]['origin'][2, 0, :]
                if fc is not None:
                     b, a = sig.butter(2, fc / (fs / 2), 'low')
                     z = sig.filtfilt(b, a, z)  
                # Differentiate
                zp = np.gradient(z)
                
                signal=np.vstack((signal,z,zp))
                

                # Using peakfind to get highest peaks in 2nd derivative of z
                mDown, _ = find_peaks(-z, distance=0.5*fs*tmpcor)
                aUp, _ = find_peaks(zp, distance=0.5*fs*tmpcor)

                # Make structured output
                if 'right' in str(traj["segment"][0,i_seg]["name"]).lower():
                    events['rto'] = aUp
                    events['rhs'] = mDown
                else:
                    events['lto'] = aUp
                    events['lhs'] = mDown
        signal=signal.T   
        signal=signal[:,1:5]                    
    else:
        if fc is not None:
             b, a = sig.butter(2, fc / (fs / 2), 'low')
             signal = sig.filtfilt(b, a, signal.T).T  
        # Get some temporary heelstrikes, to get an idea of main frequency
        nEst = int(min(120 * fs, len(signal[:, 0])))

        # Check power spectrum for first peak
        f, pxx = sig.welch(signal[:nEst, 1], fs, nperseg=2048, noverlap=25, nfft=2048, scaling='density')
        try:
            M = np.argmax(pxx)
            tmp = (f[M] * fs)*1.7
        except:
            rhs = np.where(np.diff(signal[:nEst, 1] > 0) == 1)[0]#TODO
            tmp = np.mean(np.diff(rhs))

        # If specified use the template correction
        tmp *= tmpcor
      
        # Find heel contacts and toe contacts, using main period
        maxout,tmp1  = find_peaks(signal[:, 1], distance=tmp)
        minout,tmp1 = find_peaks(-signal[:, 1], distance=tmp)
        hc = minout
        to = maxout

        # Create y-signal with no drift
        b, a = sig.butter(4, 0.5 / (fs / 2), 'low')
        y = sig.filtfilt(b, a, signal[:, 0])
        y = signal[:, 0] - y

        # Split heel and toe to left and right
        events['lhs'] = hc[y[hc] > 0].astype(int)
        events['rhs'] = hc[y[hc] < 0].astype(int)
        events['lto'] = to[y[to] > 0].astype(int)
        events['rto'] = to[y[to] < 0].astype(int)
    try:
        events,flag = order_events(events)
    except:
        print('did not order events for some reason')

    events['fs'] = fs
    return events, signal


def spline_interp_find_gaps(z, gap_size):
    # TODO This function needs to be implemented
    pass




def order_events(events, loc_mode='walk'):
    """
    Order gait events and check for consistency.

    Parameters:
    events (dict): Dictionary containing event indices
    loc_mode (str): 'walk' or 'run' mode

    Returns:
    tuple: (events, flag)
        events (dict): Ordered event indices
        flag (int): 0 if ordered correctly, 1 or 2 if issues detected
    """

    lhs = events['lhs']
    rhs = events['rhs']
    lto = events['lto']
    rto = events['rto']

    if loc_mode == 'walk':
        event_list = np.concatenate([lhs, rto, rhs, lto])
        lhs_code = np.ones(len(lhs)) * 1
        rto_code = np.ones(len(rto)) * 2
        rhs_code = np.ones(len(rhs)) * 3
        lto_code = np.ones(len(lto)) * 4
        code = np.concatenate([lhs_code, rto_code, rhs_code, lto_code])

        # Make one data matrix, sort
        data = np.column_stack((code, event_list))
        data = data[data[:, 1].argsort()]

        # Throw away everything before first left hs, and after last lto
        ind_begin = np.where(data[:, 0] == 1)[0][0]
        ind_end = np.where(data[:, 0] == 4)[0][-1]
        data = data[ind_begin:ind_end+1, :]

        # Rewrite the lhs, lto, etc.
        events['lhs'] = data[data[:, 0] == 1, 1].astype(int)
        events['rto'] = data[data[:, 0] == 2, 1].astype(int)
        events['rhs'] = data[data[:, 0] == 3, 1].astype(int)
        events['lto'] = data[data[:, 0] == 4, 1].astype(int)

        # Check if same amount of each event, and if order is correct
        flag = 0
        if not (len(events['lto']) == len(events['rto']) == len(events['rhs']) == len(events['lhs'])):
            flag = 1
        elif (any(events['rto'] - events['lhs'] < 0) or
              any(events['rhs'] - events['rto'] < 0) or
              any(events['lto'] - events['rhs'] < 0) or
              any(events['lhs'][1:] - events['lto'][:-1] < 0)):
            flag = 2

    elif loc_mode == 'run':
        event_list = np.concatenate([lhs, lto, rhs, rto])
        lhs_code = np.ones(len(lhs)) * 1
        lto_code = np.ones(len(lto)) * 2
        rhs_code = np.ones(len(rhs)) * 3
        rto_code = np.ones(len(rto)) * 4
        code = np.concatenate([lhs_code, lto_code, rhs_code, rto_code])

        # Make one data matrix, sort
        data = np.column_stack((code, event_list))
        data = data[data[:, 1].argsort()]

        # Throw away everything before first left hs, and after last lto
        ind_begin = np.where(data[:, 0] == 1)[0][0]
        ind_end = np.where(data[:, 0] == 4)[0][-1]
        data = data[ind_begin:ind_end+1, :]

        # Rewrite the lhs, lto, etc.
        events['lhs'] = data[data[:, 0] == 1, 1].astype(int)
        events['lto'] = data[data[:, 0] == 2, 1].astype(int)
        events['rhs'] = data[data[:, 0] == 3, 1].astype(int)
        events['rto'] = data[data[:, 0] == 4, 1].astype(int)

        flag = 0
        if not (len(events['lto']) == len(events['rto']) == len(events['rhs']) == len(events['lhs'])):
            flag = 1

    return events, flag




def normalizetimebase(signal, trigind=None):
    """
    Normalizes timebase and ensemble-averages cyclic signals

    Parameters:
    signal (numpy.array): Any one-dimensional array
    trigind (numpy.array): Array of indices, for instance heelstrike moments

    Returns:
    tuple: (Cycle, TimeGain)
        Cycle (numpy.array): Normalized and ensemble-averaged signal
        TimeGain (float): (average) amplification/reduction of time-axis
    """

    if trigind is None:
        trigind = np.array([0, len(signal) - 1])

    nansignal = np.zeros_like(signal)
    nansignal[~np.isnan(signal)] = 1

    Cycle = np.full((101, 1), np.nan)
    Cyclenan = np.full((101, 1), np.nan)
    # CycleLength = -101

    N = len(trigind) - 1

    if N > 1:
        Cycle = np.full((101, N + 2), np.nan)
        Cyclenan = np.full((101, N), np.nan)
        CycleLength = np.zeros(N)

        for i in range(N):
            x = np.arange(trigind[i], trigind[i + 1]) - trigind[i]
            CycleLength[i] = len(x)
            x = x * 101 / (trigind[i + 1] - trigind[i])
            x = x + 1

            # try:
            f = interpolate.interp1d(x, signal[0,   trigind[i]:trigind[i + 1]], kind='cubic')
            Cycle[:, i] = f(np.arange(1, 102))

            f_nan = interpolate.interp1d(x, nansignal[0,trigind[i]:trigind[i + 1]], kind='cubic')
            Cyclenan[:, i] = f_nan(np.arange(1, 102))
            # except:
            #     # Cyclenan[:, i] = np.nan
            #     # Cycle[:, i] = np.nan

        Cyclenan[Cyclenan < 1] = np.nan
        # Cycle = Cycle * Cyclenan  # Uncomment if needed

        Cycle[:, N] = np.nanmean(Cycle[:, :N], axis=1)
        Cycle[:, N + 1] = np.nanstd(Cycle[:, :N], axis=1)

        TimeGain = 101 / np.mean(CycleLength)

    elif N == 1:
        x = np.arange(trigind[0], trigind[1]) - trigind[0]
        CycleLength = len(x)
        x = x * 100 / (trigind[1] - trigind[0]) + 1

        try:
            f = interpolate.interp1d(x, signal[trigind[0]:trigind[1]], kind='cubic')
            Cycle[:, 0] = f(np.arange(1, 102))
        except:
            Cycle[:, 0] = np.nan

        TimeGain = 101 / CycleLength

    return Cycle, TimeGain


    