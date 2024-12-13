# Turn off vertical migration for now; start with just physics

# Same as V1, except adding deactivation after 90 days.  See the update method

#import datetime
#import numpy as np
#import logging; logger = logging.getLogger(__name__)
#from opendrift.models.oceandrift import Lagrangian3DArray, OceanDrift
#from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED

import sys
from datetime import datetime, timedelta
import numpy as np
#from scipy.interpolate import interp1d
import logging; logger = logging.getLogger(__name__)
from opendrift.models.oceandrift import Lagrangian3DArray, OceanDrift
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED
from scipy import spatial



#---------------------------------------------------------------------------------
# Swim data
#---------------------------------------------------------------------------------
swim_data_file = '/home/blaughli/tracking_project/practice/onshore_swim_work/z_production/swim_data_simple_v2.npz'

d = np.load(swim_data_file)
onshore_swim_component_x = d['onshore_swim_component_x']
onshore_swim_component_y = d['onshore_swim_component_y']
coord_array = d['coord_array']

model_lat_avg = 38 # CHANGE IF WE NEED TO BE MORE ACCURATE (ie make specific to each particle... ugh!)

meters_per_degree_lat = 111111
meters_per_degree_lon = meters_per_degree_lat * np.cos(np.radians(model_lat_avg))

# Also for now, just set a constant max swim speed (later should make it dependent on larvae size)
# From Kashef paper, speeds for rockfish larvae: at parturation: 0.5-1.8 cm/s;"newly settled": 8.6 to 53.5cm/s

#swim_speed_max = 0 # 0 cm/s  ???? Why do i keep stranding particles?
#swim_speed_max = 0.0001 #.1cm/s  ???? Why do i keep stranding particles?
#swim_speed_max = 0.001 #.1cm/s  ???? Why do i keep stranding particles?
#swim_speed_max = 0.01 #1cm/s  ???? Why do i keep stranding particles?
swim_speed_max = 0.1 #10cm/s
#swim_speed_max = 10 #10cm/s  WHY NOT, nothing is working anyway




#---------------------------------------------------------------------------------
# At some point may want to define a class for the particles, like:

# class LarvalElement(Lagrangian3DArray):
#---------------------------------------------------------------------------------

# Define the particle deactivation time (seconds).
#drift_days = 30
#drift_days = 10
#drift_days = 5
#drift_days = 150
drift_days = 180
drift_seconds = (drift_days) * 24 * 60 * 60 


class LarvalDispersal(OceanDrift):
    """Following example of LarvalFish, and trying to add behavior for larva

    """
    # Again, may want to implement and then use your own class of element, a-la:
    #ElementType = LarvalElement
    ElementType = Lagrangian3DArray

    #max_speed = 1  # m/s     # Why was this here??? Did I have a specific reason, or did I copy it from OceanDrift???

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_height': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'upward_sea_water_velocity': {'fallback': 0},
        'ocean_vertical_diffusivity': {'fallback': 0.00001, 'profiles': True},
        'surface_downward_x_stress': {'fallback': 0},
        'surface_downward_y_stress': {'fallback': 0},
        'turbulent_kinetic_energy': {'fallback': 0},
        'turbulent_generic_length_scale': {'fallback': 0},
        'ocean_mixed_layer_thickness': {'fallback': 50},
        'sea_floor_depth_below_sea_level': {'fallback': 10000},
        'land_binary_mask': {'fallback': None},
        'sea_water_temperature': {'fallback': 10, 'profiles': True},
        'sea_water_salinity': {'fallback': 34, 'profiles': True},
        'CalC': {'fallback': 0},
        'DON': {'fallback': 0},
        'NH4': {'fallback': 0},
        'NO3': {'fallback': 0},
        'PON': {'fallback': 0},
        'Pzooplankton': {'fallback': 0},
        'SiOH4': {'fallback': 0},
        'TIC': {'fallback': 0},
        'alkalinity': {'fallback': 0},
        'diatom': {'fallback': 0},
        'mesozooplankton': {'fallback': 0},
        'microzooplankton': {'fallback': 0},
        'nanophytoplankton': {'fallback': 0},
        'omega': {'fallback': 0},
        'opal': {'fallback': 0},
        'oxygen': {'fallback': 0},
        'pCO2': {'fallback': 0},
        'pH': {'fallback': 0},
        }

        # Deleted all having to do with waves, and not sure I need these either:
        # (see old versions for the exhaustive list)



    #required_profiles_z_range = [0, -50]  # The depth range (in m) which profiles should cover

    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(LarvalDispersal, self).__init__(*args, **kwargs)
        
        self._add_config({
            'IBM:fraction_of_timestep_swimming':
                {'type': 'float', 'default': 0.15,
                 'min': 0.0, 'max': 1.0, 'units': 'fraction',
                 'description': 'Fraction of timestep swimming',
                 'level': CONFIG_LEVEL_ADVANCED},
            'drift:max_lifespan_days':{
                    'type': 'int',
                    'default': 180,
                    'min': 1,
                    'max': 9999,
                    'units': 'days',
                    'description': 'Maximum drifter lifespan before deactivation',
                    'level': CONFIG_LEVEL_BASIC},
        })
       
        
        #self._set_config_default('behavior:onshore_swimming', False)
        ##self._set_config_default('behavior:onshore_swimming', True)

        #self._set_config_default('drift:vertical_advection', True)
        #self._set_config_default('drift:vertical_mixing', True)
        #self._set_config_default('general:coastline_action', 'previous')

        #self._set_config_default('drift:profile_depth', 50)  # The depth range (in m) which profiles should cover
        ###self._set_config_default('drift:profile_depth', [0, -50])  # The depth range (in m) which profiles should cover

    # ---------------------------------------------------------------------------------------------
    # Will we want to update properties of the larvae?  See LarvalFish for what was here

    #def update_terminal_velocity(self, Tprofiles=None,
    #                             Sprofiles=None, z_index=None):
    #    """Calculate terminal velocity for Pelagic Egg

    #    according to
    #    S. Sundby (1983): A one-dimensional model for the vertical
    #    distribution of pelagic fish eggs in the mixed layer
    #    Deep Sea Research (30) pp. 645-661

    #    Method copied from ibm.f90 module of LADIM:
    #    Vikebo, F., S. Sundby, B. Aadlandsvik and O. Otteraa (2007),
    #    Fish. Oceanogr. (16) pp. 216-228
    #    """

    
    #def update_larvae(self):
    # ---------------------------------------------------------------------------------------------

    # Copied from LarvalFish 
    def larvae_onshore_swim(self):
        
        lons = self.elements.lon.reshape(-1,1)
        lats = self.elements.lat.reshape(-1,1)
        lonlatpairs = np.concatenate((lons,lats), axis=1)
        distances,indices = spatial.KDTree(coord_array).query(lonlatpairs)

#        print('USER PRINT STATEMENT: vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv',flush=True)
#        print('USER PRINT STATEMENT: number of floats FOUND BY KDTree: {} '.format(len(indices)),flush=True)
#        print('USER PRINT STATEMENT: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',flush=True)

        print(f'index: {indices[0]}')

        self.elements.lon = self.elements.lon + onshore_swim_component_x[indices]*swim_speed_max/meters_per_degree_lon*self.time_step.total_seconds()
        self.elements.lat = self.elements.lat + onshore_swim_component_y[indices]*swim_speed_max/meters_per_degree_lat*self.time_step.total_seconds()

        # TEST
        #self.elements.lon = self.elements.lon - onshore_swim_component_x[indices]*swim_speed_max
        #self.elements.lat = self.elements.lat - onshore_swim_component_y[indices]*swim_speed_max
        
        
        # TEST
        #self.elements.lon = self.elements.lon - swim_speed_max
        #self.elements.lat = self.elements.lat - swim_speed_max

#        vx = np.ones_like(self.elements.lon)*swim_speed_max
#        vy = np.ones_like(self.elements.lon)*swim_speed_max

#        self.update_positions(vx, vy)

    def update(self):
        """Update positions and properties of elements."""
        # copied from "OceanDrift", with the addition of the deactivation

        # Simply move particles with ambient current
#        self.advect_ocean_current()

        # Advect particles due to surface wind drag,
        # according to element property wind_drift_factor
        #self.advect_wind()

        # Stokes drift
        #self.stokes_drift()

#        # Turbulent Mixing
#        if self.get_config('drift:vertical_mixing') is True:
#            self.update_terminal_velocity()
#            self.vertical_mixing()
#        else:  # Buoyancy
#            self.vertical_buoyancy()

#        # Vertical advection
#        if self.get_config('drift:vertical_advection') is True:
#            self.vertical_advection()
        
        ### Swimming
        ##self.get_config('behavior:onshore_swimming' is True:
        self.larvae_onshore_swim()


        #vx = np.ones_like(self.elements.lon)*swim_speed_max
        #vy = np.ones_like(self.elements.lon)*swim_speed_max

        #vx = np.ones(1)*10
        #vy = np.ones(1)*10

        #self.update_positions(vx, vy)
        
        # TEST
        #self.elements.lon = self.elements.lon + 0.1
        #self.elements.lat = self.elements.lat + 0.1


        # Deactivate floats after "drift_days" has passed
#        self.deactivate_elements(self.elements.age_seconds > drift_seconds, reason='age > {} days'.format(drift_days))



