import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import plotly.graph_objects as go

class AcousticsDataset:
    """
    A custom class to relate the outputs from different JPP generated (ROOT) files to be visualised in python.
    INPUTS:
        -> detectorfile
        -> hitfile
        -> debug (optional)
        -> triggerfile
        -> displayfile

    In addition, the constructer loads the detector geometry through running JPrintDetector on the detector file which is located in:
    path/to/detector/<detector_name>.datx

    If this file does not exist OR the environment is not suited for running JPP software, the constructor will look for a matching .txt file locally
    """

    def __init__(self, *args):
        self.d = 1

        if len(args) < 2:
            raise ValueError("Please provide 2, 3 or 5 arguments")
        else:
            detectorfile = args[0]
            hitfile      = args[1]

            self.detector_file  = detectorfile
            self.detector_name  = detectorfile.split('.')[0].split('/')[-1]

            # Load detector geometry
            self.hydrophone_dict  = self._load_detector_geometry(self.detector_file)

            # load ROOT hit files and content
            self.hitfile            = ROOT.TFile(hitfile)
            self.event_tree         = self.hitfile.Get("ACOUMCEVENT")
            self.toa_tree           = self.hitfile.Get("ACOUHIT")
            self.numevents          = self.hitfile.Get("ACOUMCEVENT").GetEntries()

            if len(args) >= 3:
                self.d = args[2]

            if len(args) == 5:
                trigfile    = args[3]
                displfile   = args[4]


                self.triggerfile      = ROOT.TFile(trigfile)
                self.displayfile      = ROOT.TFile(displfile)

                # Load trees
                self.detected_tree      = self.displayfile.Get("tree")
                self.isprocessed        = 1
            
            if (len(args) == 4 or len(args) > 5):
                raise ValueError("Please provide 2, 3 or 5 arguments")


    def _load_detector_geometry(self, detector_file):
        """Reads the DOM geometry using JPrintDetector if possible, otherwise tries local file."""
        hydrophone_dict = {}
        detector_name   = self.detector_name
        export_dir      = os.path.join(".", detector_name)
        os.makedirs(export_dir, exist_ok=True)
        export_path     = os.path.join(export_dir, "hydrophone_dict_export.txt")

        # First, try to run JPrintDetector
        try:
            self._DEBUG(f"Trying to extract DOM geometry using JPrintDetector from: {detector_file}")
            result  = subprocess.run(['JPrintDetector', '-a', detector_file], stdout=subprocess.PIPE, text=True, check=True)
            lines   = result.stdout.strip().split('\n')
            for i, line in enumerate(lines):
                if i > 9 and line.strip() != "" and not line.startswith("#"):
                    parts = line.split()
                    if len(parts) >= 6:
                        hydrophoneid    = int(parts[0])
                        x, y, z         = float(parts[3]), float(parts[4]), float(parts[5])
                        hydrophone_dict[hydrophoneid] = [x, y, z]
            # Export hydrophone_dict to a local file for reference
            with open(export_path, "w") as fout:
                for hydrophoneid, coords in hydrophone_dict.items():
                    fout.write(f"{hydrophoneid} {coords[0]} {coords[1]} {coords[2]}\n")
            self._DEBUG(f"DOM geometry successfully extracted and saved to: {export_path}")
            return hydrophone_dict
        except Exception as e:
            self._DEBUG(f"JPrintDetector failed: {e}")

        # If JPrintDetector fails, try to load from local file
        if os.path.exists(export_path):
            self._DEBUG(f"Loading DOM geometry from local export: {export_path}")
            with open(export_path, "r") as fin:
                for line in fin:
                    parts = line.strip().split()
                    if len(parts) == 4:
                        hydrophoneid = int(parts[0])
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        hydrophone_dict[hydrophoneid] = [x, y, z]
            self._DEBUG(f"DOM geometry successfully loaded from: {export_path}")
            return hydrophone_dict

        # If both methods fail, raise an error
        raise RuntimeError(f"Could not load DOM geometry using JPrintDetector or from local file: {export_path}")


    def correlate_events(self):
        if not self.isprocessed:
            raise ValueError("Trigger files not provided")
        
        """Correlates generated events with detected events based on ToA. ChatGPT optimized."""

        print("Correlating generated events with detected events")
        
        detected_eventdict = {}
        for entry in self.detected_tree:
            eid = int(entry.EventID) - 1

            did = int(entry.HydrophoneID)
            tid = int(entry.TemplateID)
            toa = entry.ToA
            snr = entry.SNR

            if eid not in detected_eventdict:
                detected_eventdict[eid] = [[], [], [], []]
            detected_eventdict[eid][0].append(did)
            detected_eventdict[eid][1].append(toa)
            detected_eventdict[eid][2].append(snr)
            detected_eventdict[eid][3].append(tid)
        
        self.detected_eventdict = detected_eventdict

        # Extract first ToA (converted to us)
        t0s = np.array([int(detected_eventdict[k][1][0] * 1e6) for k in detected_eventdict]) # TODO CHECK

        # Build generated TOA_US as a NumPy array
        T0GENUS = []
        for evt in self.event_tree:
            T0GENUS.append(evt.t0_s * 1e6)
        T0GENUS = np.array(T0GENUS)

        # Sort TOAUS for fast lookup
        sorted_indices = np.argsort(T0GENUS)
        sorted_T0GENUS = T0GENUS[sorted_indices]

        # Efficient nearest neighbor search
        self.detected_event_generation_IDs = np.zeros_like(t0s, dtype=int)
        for i, t0 in enumerate(t0s):
            idx = np.searchsorted(sorted_T0GENUS, t0)

            # Clamp to valid range
            idxs = [max(0, idx - 1), min(len(sorted_T0GENUS) - 1, idx)]
            
            # Choose the closer one
            nearest_idx = idxs[np.argmin([abs(sorted_T0GENUS[j] - t0) for j in idxs])]
            TOA_ID = int(sorted_indices[nearest_idx])

            # self.event_tree.GetEntry(TOA_ID)
            self.detected_event_generation_IDs[i] = TOA_ID # CHECK
        
        print("Correlation complete.")
        
    
    def print_tree_info(self):
        """Prints the contents of the ROOT file and branches of all TTrees."""

        print("Keys in the ROOT file:")
        for key in self.hitfile.GetListOfKeys():
            print(f"  - {key.GetName()} ({key.GetClassName()})")

        print("\nListing branches of TTree objects:")
        for key in self.hitfile.GetListOfKeys():
            obj = self.hitfile.Get(key.GetName())
            if obj.InheritsFrom("TTree"):
                print(f"\nTree: {key.GetName()}")
                obj.Print()


    def plot_TOA_histogram(self):
        """Plots a histogram of Arrival timnes from the TOA tree."""
        toa_us = [entry.TOA_US for entry in self.toa_tree]
        plt.figure(figsize=(10, 6))
        plt.hist(toa_us, bins=1000, histtype='step', color='green')
        plt.xlabel(r"Toa [$\mu$s]")
        plt.ylabel("Counts")
        plt.title("Histogram of arrival times")
        plt.grid(True)
        plt.show()


    def plot_all_events(self):
        """Plots all events with direction arrows."""
        Nentries = self.event_tree.GetEntries()
        x, y, z, dx, dy, dz = np.zeros(Nentries), np.zeros(Nentries), np.zeros(Nentries), np.zeros(Nentries), np.zeros(Nentries), np.zeros(Nentries)

        for i, entry in enumerate(self.event_tree):
            x[i], y[i], z[i]    = entry.x,  entry.y,  entry.z
            dx[i], dy[i], dz[i] = entry.dx, entry.dy, entry.dz

        fig   = plt.figure(figsize=(10, 8))
        ax    = fig.add_subplot(111, projection='3d')
        ax    .quiver(x, y, z, dx, dy, dz, length=150, normalize=True, color='k', linewidth=1)
        sc    = ax.scatter(x, y, z, c=dz, cmap='coolwarm', s=20)
        cbar  = plt.colorbar(sc, ax=ax, pad=0.1)
        cbar  .set_label('dz')
        ax    .set_xlabel('x')
        ax    .set_ylabel('y')
        ax    .set_zlabel('z')
        ax    .set_title('3D Scatter Plot with Direction Arrows')
        plt   .show()


    def plot_event(self, eventID=None):
        """Plots the 3D hydrophone hits for a single event."""
        if eventID is None:
            eventID = np.random.randint(0, self.numevents) 
        
        if eventID < 0 or eventID >= self.numevents:
            raise ValueError(f"Invalid event ID: {eventID}. Must be between 0 and {self.numevents - 1}.")
        
        self.event_tree.GetEntry(eventID)
        x, y, z     = self.event_tree.x,  self.event_tree.y,  self.event_tree.z
        dx, dy, dz  = self.event_tree.dx, self.event_tree.dy, self.event_tree.dz
        TGEN_US     = self.event_tree.t0_s * 1e6 
        print(rf"Plotting event ID {eventID} at position ({x:.0f}, {y:.0f}, {z:.0f}) with direction ({dx:.3f}, {dy:.3f}, {dz:.3f}) at t = {TGEN_US:.0f} $\mu$s")
        
        # Get hits for this event from TOA tree inside a 5 sec window
        HydrophoneIDs, TOA_US = [], []
        for entry in self.toa_tree:
            if (entry.TOA_US > TGEN_US and entry.TOA_US < TGEN_US + 5e6):
                HydrophoneIDs.append(int(entry.HydrophoneID))
                TOA_US.append(entry.TOA_US)
        if not HydrophoneIDs:
            raise ValueError(f"No hits for event ID {eventID}")
    
        all_HydrophoneIDs = np.array(list(self.hydrophone_dict.keys()))
        all_positions     = np.array([self.hydrophone_dict[hydrophoneid] for hydrophoneid in all_HydrophoneIDs])

        hit_HydrophoneIDs = np.array(HydrophoneIDs)
        hit_positions     = np.array([self.hydrophone_dict[hydrophoneid] for hydrophoneid in hit_HydrophoneIDs])
        hit_toa_us        = np.array(TOA_US) / 1e6

        nonhit_mask = ~np.isin(all_HydrophoneIDs, hit_HydrophoneIDs)
        nonhit_positions = all_positions[nonhit_mask]

        fig = go.Figure()
        fig.add_trace(go.Scatter3d(x=nonhit_positions[:, 0], y=nonhit_positions[:, 1], z=nonhit_positions[:, 2],
                                   mode='markers', marker=dict(size=2, color='black', opacity=0.15)))

        fig.add_trace(go.Scatter3d(x=hit_positions[:, 0], y=hit_positions[:, 1], z=hit_positions[:, 2],
                                   mode='markers', marker=dict(size=5, color=hit_toa_us, colorscale='Plasma',
                                                               colorbar=dict(title='Toa [s]'))))

        fig.add_trace(go.Scatter3d(x=[x], y=[y], z=[z], mode='markers',
                                   marker=dict(size=10, color='lime', symbol='cross')))

        direction     = np.array([dx, dy, dz])
        direction     = direction / np.linalg.norm(direction) if np.linalg.norm(direction) else np.array([0, 0, 1])
        shaft_length  = 1000
        head_length   = 1000
        start_point   = np.array([x, y, z])
        end_point     = start_point + shaft_length * direction
        cone_base     = end_point
        cone_tip      = cone_base + head_length * direction

        fig.add_trace(go.Scatter3d(x=[start_point[0], end_point[0]], y=[start_point[1], end_point[1]],
                                   z=[start_point[2], end_point[2]], mode='lines',
                                   line=dict(color='lime', width=6)))

        fig.add_trace(go.Cone(x=[cone_base[0]], y=[cone_base[1]], z=[cone_base[2]],
                              u=[direction[0]], v=[direction[1]], w=[direction[2]],
                              sizemode='absolute', sizeref=500, anchor='tail', showscale=False,
                              colorscale=[[0, 'lime'], [1, 'lime']]))

        fig.update_layout(width=1000, height=800,
                          margin=dict(l=0, r=0, t=50, b=0),
                          title=f'3D Hydrophone Hits for Event ID {eventID}',
                          showlegend=False)
        fig.show()


    def plot_detected_event(self, detected_eventID=None):
        """Plots the 3D hydrophone hits for a detected event."""

        if not self.isprocessed:
            raise ValueError("Trigger files not provided")
        
        # Relate the events from the detection to the events from the generation
        if not hasattr(self, 'detected_event_generation_IDs'):
            self.correlate_events()
        
        if detected_eventID is None:
            maxID = int(np.unique(self.displayfile.Get("tree").EventID)[0])
            detected_eventID = np.random.randint(0, maxID) 
        
        generated_eventID = int(self.detected_event_generation_IDs[detected_eventID])

        self.event_tree.GetEntry(generated_eventID)
        x, y, z     = self.event_tree.x,  self.event_tree.y,  self.event_tree.z
        dx, dy, dz  = self.event_tree.dx, self.event_tree.dy, self.event_tree.dz
        TGEN_US     = self.event_tree.t0_s * 1e6 

        print(f"Plotting event ID {generated_eventID} at position ({x:.0f}, {y:.0f}, {z:.0f}) with direction ({dx:.3f}, {dy:.3f}, {dz:.3f})")

        # Get hits for this event from TOA tree
        HydrophoneIDs, TOA_US = [], []
        for entry in self.toa_tree:
            if (entry.TOA_US > TGEN_US and entry.TOA_US < TGEN_US + 8e6):
                HydrophoneIDs.append(entry.HydrophoneID)
                TOA_US.append(entry.TOA_US)
        if not HydrophoneIDs:
            raise ValueError(f"No hits for event ID {detected_eventID}")

        # Get the HydrophoneIDs and hits from DETECTION

        all_HydrophoneIDs       = np.array(list(self.hydrophone_dict.keys()))
        all_positions           = np.array([self.hydrophone_dict[hydrophoneid] for hydrophoneid in all_HydrophoneIDs])

        hit_HydrophoneIDs       = np.array(HydrophoneIDs)
        hit_positions           = np.array([self.hydrophone_dict[hydrophoneid] for hydrophoneid in hit_HydrophoneIDs])
        hit_toa_us              = np.array(TOA_US) / 1e9
        
        detected_HydrophoneIDs  = np.array(self.detected_eventdict[detected_eventID][0])
        detected_hit_positions  = np.array([self.hydrophone_dict[hydrophoneid] for hydrophoneid in detected_HydrophoneIDs])

        nonhit_mask             = ~np.isin(all_HydrophoneIDs, hit_HydrophoneIDs)
        nonhit_positions        = all_positions[nonhit_mask]

        # Print the # of doms used in detection
        print(f"{len(detected_HydrophoneIDs)} out of {len(hit_HydrophoneIDs)} used in JBillabong: {len(detected_HydrophoneIDs)/len(hit_HydrophoneIDs)*100:.1f}%")
        
        fig = go.Figure()
        fig.add_trace(go.Scatter3d(x=nonhit_positions[:, 0], y=nonhit_positions[:, 1], z=nonhit_positions[:, 2],
                                    mode='markers', marker=dict(size=2, color='black', opacity=0.15)))

        fig.add_trace(go.Scatter3d(x=hit_positions[:, 0], y=hit_positions[:, 1], z=hit_positions[:, 2],
                                    mode='markers', marker=dict(size=5, color=hit_toa_us, colorscale='Plasma',
                                                                colorbar=dict(title='Toa [s]'))))

        fig.add_trace(go.Scatter3d(x=detected_hit_positions[:, 0], y=detected_hit_positions[:, 1], z=detected_hit_positions[:, 2],
                                    mode='markers', marker=dict(size=10, color='black', symbol='cross')))

        fig.add_trace(go.Scatter3d(x=[x], y=[y], z=[z], mode='markers',
                                    marker=dict(size=10, color='lime', symbol='cross')))

        direction     = np.array([dx, dy, dz])
        direction     = direction / np.linalg.norm(direction) if np.linalg.norm(direction) else np.array([0, 0, 1])
        shaft_length  = 1000
        head_length   = 1000
        start_point   = np.array([x, y, z])
        end_point     = start_point + shaft_length * direction
        cone_base     = end_point
        cone_tip      = cone_base + head_length * direction

        fig.add_trace(go.Scatter3d(x=[start_point[0], end_point[0]], y=[start_point[1], end_point[1]],
                                    z=[start_point[2], end_point[2]], mode='lines',
                                    line=dict(color='lime', width=6)))

        fig.add_trace(go.Cone(x=[cone_base[0]], y=[cone_base[1]], z=[cone_base[2]],
                                u=[direction[0]], v=[direction[1]], w=[direction[2]],
                                sizemode='absolute', sizeref=500, anchor='tail', showscale=False,
                                colorscale=[[0, 'lime'], [1, 'lime']]))

        fig.update_layout(width=1000, height=800,
                            margin=dict(l=0, r=0, t=50, b=0),
                            title=f'3D Hydrophone Hits for generated Event ID {generated_eventID}, Detected Event ID {detected_eventID}',
                            showlegend=False)
        fig.show()

    def compute_detected_directions(self):
        """computes angular distribution of detected events."""

        if not self.isprocessed:
            raise ValueError("Trigger files not provided")
        
        # Relate the events from the detection to the events from the generation
        if not hasattr(self, 'detected_event_generation_IDs'):
            self.correlate_events()
        
        dzs = np.zeros(len(self.detected_event_generation_IDs))
        for i, generated_eventID in enumerate(self.detected_event_generation_IDs):
            self.event_tree.GetEntry(int(generated_eventID))
            
            dz = self.event_tree.dz
            dzs[i]  = dz
        
        return dzs
        
    def compute_generated_directions(self):
        """computes angular distribution of generated events."""
        dzs = np.zeros(self.numevents)

        for i in range(self.numevents):
            self.event_tree.GetEntry(int(i))
            
            dz = self.event_tree.dz
            dzs[i]  = dz
        
        return dzs

    def _what_hydrophones(self):
        """Compute what hydrophones are used and how often"""
        print("Computing hydrophones-used histogram")
        used_hydro_hist = {}

        # Loop over all events in self.detected_eventdict
        for i in self.detected_eventdict:
            hydro_ids = self.detected_eventdict[i][0]
            
            for hid in hydro_ids:
                # Increase count or initialize to 1
                used_hydro_hist[hid] = used_hydro_hist.get(hid, 0) + 1
        
        self.used_hydro_hist = used_hydro_hist
 

    def plot_what_hydrophones_are_used(self):
        """computes what hydrophones are actually being used."""
        if not self.isprocessed:
            raise ValueError("Trigger files not provided")
        
        if not hasattr(self, 'used_hydro_hist'):

            # Relate the events from the detection to the events from the generation
            if not hasattr(self, 'detected_event_generation_IDs'):
                self.correlate_events()
            
            self._what_hydrophones()
        
        # Now for the plotting
        all_HydrophoneIDs       = np.array(list(self.hydrophone_dict.keys()))
        all_positions           = np.array([self.hydrophone_dict[hydrophoneid] for hydrophoneid in all_HydrophoneIDs])

        used_HydrophoneIDs      = np.array(list(self.used_hydro_hist.keys()))
        used_positions          = np.array([self.hydrophone_dict[hydrophoneid] for hydrophoneid in used_HydrophoneIDs])
        frequency               = np.array([self.used_hydro_hist[hydrophoneid] for hydrophoneid in used_HydrophoneIDs])
        
        nonhit_mask             = ~np.isin(all_HydrophoneIDs, used_HydrophoneIDs)
        nonhit_positions        = all_positions[nonhit_mask]

        fig = go.Figure()
        fig.add_trace(go.Scatter3d(x=nonhit_positions[:, 0], y=nonhit_positions[:, 1], z=nonhit_positions[:, 2],
                                    mode='markers', marker=dict(size=2, color='black', opacity=0.15)))

        fig.add_trace(go.Scatter3d(x=used_positions[:, 0], y=used_positions[:, 1], z=used_positions[:, 2],
                                    mode='markers', marker=dict(size=5, color=frequency, colorscale='Plasma',
                                                                colorbar=dict(title='Frequency'))))
        
        fig.update_layout(width=1000, height=800,
                            margin=dict(l=0, r=0, t=50, b=0),
                            title=f'What hydrophones are used',
                            showlegend=False)
        fig.show()           


    def all_detected_vertices(self):
        """Compute the detected event vertices"""
        if not self.isprocessed:
            raise ValueError("Trigger files not provided")

        # Relate the events from the detection to the events from the generation
        if not hasattr(self, 'detected_event_generation_IDs'):
            self.correlate_events()

        if len(self.detected_event_generation_IDs) == 0:
            return np.array([]), np.array([]), np.array([])
        
        X = np.zeros_like(self.detected_event_generation_IDs)
        Y = np.zeros_like(self.detected_event_generation_IDs)
        Z = np.zeros_like(self.detected_event_generation_IDs)

        for i, generated_eventID in enumerate(self.detected_event_generation_IDs):

            self.event_tree.GetEntry(int(generated_eventID))
            X[i] = self.event_tree.x
            Y[i] = self.event_tree.y
            Z[i] = self.event_tree.z
        
        return X,Y,Z


    def _DEBUG(self,txt):
        if self.d: print(txt)


######################################################################################################################
######################################                  NOISE                   ######################################
######################################################################################################################

class AcousticsNoiseDataset:
    """
    A custom class to relate the outputs from different JPP generated (ROOT) files to be visualised in python.
    INPUTS:
        -> detectorfile
        -> hitfile
        -> debug (optional)
        -> triggerfile
        -> displayfile

    In addition, the constructer loads the detector geometry through running JPrintDetector on the detector file which is located in:
    path/to/detector/<detector_name>.datx

    If this file does not exist OR the environment is not suited for running JPP software, the constructor will look for a matching .txt file locally
    """

    def __init__(self, *args):
        self.d = 1

        if len(args) < 2:
            raise ValueError("Please provide 2, 3 or 5 arguments")
        else:
            detectorfile = args[0]
            hitfile      = args[1]

            self.detector_file  = detectorfile
            self.detector_name  = detectorfile.split('.')[0].split('/')[-1]

            # Load detector geometry
            self.hydrophone_dict  = self._load_detector_geometry(self.detector_file)

            # load ROOT hit files and content
            self.hitfile          = ROOT.TFile(hitfile)
            self.event_tree         = self.hitfile.Get("ACOUMCEVENT")
            self.toa_tree           = self.hitfile.Get("ACOUHIT")

            if len(args) >= 3:
                self.d = args[2]

            if len(args) == 5:
                trigfile    = args[3]
                displfile   = args[4]


                self.triggerfile      = ROOT.TFile(trigfile)
                self.displayfile      = ROOT.TFile(displfile)

                # Load trees
                self.detected_tree      = self.displayfile.Get("tree")
                self.isprocessed        = 1
            
            if (len(args) == 4 or len(args) > 5):
                raise ValueError("Please provide 2, 3 or 5 arguments")


    def _load_detector_geometry(self, detector_file):
        """Reads the DOM geometry using JPrintDetector if possible, otherwise tries local file."""
        hydrophone_dict = {}
        detector_name   = self.detector_name
        export_dir      = os.path.join(".", detector_name)
        os.makedirs(export_dir, exist_ok=True)
        export_path     = os.path.join(export_dir, "hydrophone_dict_export.txt")

        # First, try to run JPrintDetector
        try:
            self._DEBUG(f"Trying to extract DOM geometry using JPrintDetector from: {detector_file}")
            result  = subprocess.run(['JPrintDetector', '-a', detector_file], stdout=subprocess.PIPE, text=True, check=True)
            lines   = result.stdout.strip().split('\n')
            for i, line in enumerate(lines):
                if i > 9 and line.strip() != "" and not line.startswith("#"):
                    parts = line.split()
                    if len(parts) >= 6:
                        hydrophoneid    = int(parts[0])
                        x, y, z         = float(parts[3]), float(parts[4]), float(parts[5])
                        hydrophone_dict[hydrophoneid] = [x, y, z]
            # Export hydrophone_dict to a local file for reference
            with open(export_path, "w") as fout:
                for hydrophoneid, coords in hydrophone_dict.items():
                    fout.write(f"{hydrophoneid} {coords[0]} {coords[1]} {coords[2]}\n")
            self._DEBUG(f"DOM geometry successfully extracted and saved to: {export_path}")
            return hydrophone_dict
        except Exception as e:
            self._DEBUG(f"JPrintDetector failed: {e}")

        # If JPrintDetector fails, try to load from local file
        if os.path.exists(export_path):
            self._DEBUG(f"Loading DOM geometry from local export: {export_path}")
            with open(export_path, "r") as fin:
                for line in fin:
                    parts = line.strip().split()
                    if len(parts) == 4:
                        hydrophoneid = int(parts[0])
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        hydrophone_dict[hydrophoneid] = [x, y, z]
            self._DEBUG(f"DOM geometry successfully loaded from: {export_path}")
            return hydrophone_dict

        # If both methods fail, raise an error
        raise RuntimeError(f"Could not load DOM geometry using JPrintDetector or from local file: {export_path}")


    def plot_detected_event(self):
        """Plots the 3D hydrophone hits for a detected event."""

        if not self.isprocessed:
            raise ValueError("Trigger files not provided")

        HydrophoneIDs, TOA_US, SNRS = [], [], []

        for entry in self.displayfile.Get("tree"):
            if entry.EventID == 1:
                HydrophoneIDs.append(entry.HydrophoneID)
                TOA_US.append(entry.ToA)
                SNRS.append(entry.SNR)

        all_HydrophoneIDs       = np.array(list(self.hydrophone_dict.keys()))
        all_positions           = np.array([self.hydrophone_dict[hydrophoneid] for hydrophoneid in all_HydrophoneIDs])

        hit_HydrophoneIDs       = np.array(HydrophoneIDs)
        hit_positions           = np.array([self.hydrophone_dict[hydrophoneid] for hydrophoneid in hit_HydrophoneIDs])
        hit_toa_us              = np.array(TOA_US) 

        nonhit_mask             = ~np.isin(all_HydrophoneIDs, hit_HydrophoneIDs)
        nonhit_positions        = all_positions[nonhit_mask]

        fig = go.Figure()
        fig.add_trace(go.Scatter3d(x=nonhit_positions[:, 0], y=nonhit_positions[:, 1], z=nonhit_positions[:, 2],
                                    mode='markers', marker=dict(size=2, color='black', opacity=0.15)))

        fig.add_trace(go.Scatter3d(x=hit_positions[:, 0], y=hit_positions[:, 1], z=hit_positions[:, 2],
                                    mode='markers', marker=dict(size=5, color=hit_toa_us, colorscale='Plasma',
                                                                colorbar=dict(title='Toa [s]'))))

        fig.update_layout(width=1000, height=800,
                            margin=dict(l=0, r=0, t=50, b=0),
                            showlegend=False)
        fig.show()

    
    def _DEBUG(self,txt):
        if self.d: print(txt)