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

# Copied from "larvalfish.py" in opendrift models
class LarvalFishElement(Lagrangian3DArray):
    """
    Extending Lagrangian3DArray with specific properties for larval and juvenile stages of fish
    """

    variables = Lagrangian3DArray.add_variables([
        ('length', {'dtype': np.float32,
                    'units': 'mm',
                    'default': 0}),
        ('weight', {'dtype': np.float32,
                    'units': 'mg',
                    'default': 0.08})])



class LarvalDispersal(OceanDrift):
    """Following example of LarvalFish, and trying to add behavior for larva

    """
    # Again, may want to implement and then use your own class of element, a-la:
    ElementType = LarvalFishElement
    #ElementType = Lagrangian3DArray

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

        ## IBM configuration options
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


        #self._set_config_default('drift:vertical_advection', True)
        #self._set_config_default('drift:vertical_mixing', True)
        #self._set_config_default('general:coastline_action', 'previous')

        #self._set_config_default('drift:profile_depth', 50)  # The depth range (in m) which profiles should cover
        ###self._set_config_default('drift:profile_depth', [0, -50])  # The depth range (in m) which profiles should cover


    def fish_growth(self, weight, temperature):
        # Weight in milligrams, temperature in celcius
        # Daily growth rate in percent according to
        # Folkvord, A. 2005. "Comparison of Size-at-Age of Larval Atlantic Cod (Gadus Morhua)
        # from Different Populations Based on Size- and Temperature-Dependent Growth
        # Models." Canadian Journal of Fisheries and Aquatic Sciences.
        # Journal Canadien Des Sciences Halieutiques # et Aquatiques 62(5): 1037-52.
        GR = 1.08 + 1.79 * temperature - 0.074 * temperature * np.log(weight) \
             - 0.0965 * temperature * np.log(weight) ** 2 \
             + 0.0112 * temperature * np.log(weight) ** 3

        # Growth rate(g) converted to milligram weight (gr_mg) per timestep:
        g = (np.log(GR / 100. + 1)) * self.time_step.total_seconds()/86400
        return weight * (np.exp(g) - 1.)


    def update_fish_larvae(self):

        # Increasing weight of larvae
        avg_weight_before = self.elements.weight.mean()
        growth = self.fish_growth(self.elements.weight,
                                  self.environment.sea_water_temperature)
        self.elements.weight += growth
        avg_weight_after = self.elements.weight.mean()
        logger.debug('Growing larvae from average size %s to %s' %
              (avg_weight_before, avg_weight_after))

        # Increasing length of larvae, according to Folkvord (2005)
        w = self.elements.weight
        self.elements.length = np.exp(2.296 + 0.277 * np.log(w) - 0.005128 *np.log10(w)**2)




    # Copied from LarvalFish 
    def larvae_vertical_migration(self):

        # Note: length is in mm:
        # https://opendrift.github.io/_modules/opendrift/models/larvalfish.html#LarvalFish.larvae_vertical_migration


        # Vertical migration of Larvae
        # Swim function from Peck et al. 2006
        # Note: is "fraction of timestep swimming" ever changed?  Or stuck at 15, or whatever we set above?

        ## INIITIALLY, let's try this without growth - just set size to 15mm (average size based on Kashef 2014)
        #avg_length = 15.  # (mm)
        #L = np.ones_like(self.elements.z) * avg_length


        # Wait... units of speed can't be m/s !!!!  tiny larvae swimming 3 meters in a second????
        ##### This function doesn't make sense - it changes sign around length = 5 (so 1mm swims at -3m/s, 5mm at 0, 10mm at 3...)
        #swim_speed = (0.261*(L**(1.552*L**(-0.08))) - 5.289/L) / 1000
        
        # So, for initial tests, just use 3m/s, since that seems like a roughly average (?) value 

        avg_speed = 0.003  #????  mm/s ? # Cursory tests show that this yields ~26m of vertical movement  in 12 hours

        swim_speed = np.ones_like(self.elements.z) * avg_speed  
        
        f = self.get_config('IBM:fraction_of_timestep_swimming')
        max_migration_per_timestep = f*swim_speed*self.time_step.total_seconds()

        # Using here UTC hours. Should be changed to local solar time,
        # although a phase shift of some hours should not make much difference
        if self.time.hour < 12:
            direction = -1  # Swimming down when light is increasing
        else:
            direction = 1  # Swimming up when light is decreasing

        #print(f'Current hour: {self.time.hour}')
        #print(f'Current z: {self.elements.z[0]}')
        #print(f'direction = {direction}')
        ##logger.info(f'Current z: {self.elements.z[0]}')
        ##logger.debug(f'Current z: {self.elements.z[0]}')
        self.elements.z = np.minimum(0, self.elements.z + direction*max_migration_per_timestep)
        #print(f'z adjustment: {direction*max_migration_per_timestep}')

    def update(self):
        """Update positions and properties of elements."""
        # copied from "OceanDrift", with the addition of the deactivation

        # Update the larvae size with the copied growth model - need to verify that this is working!  Do weight and length get saved,
        # and do they compare well with literature for west coast rock fish (model is for Atlantic cod)?
        ##### FOR NOW, don't do growth.  Just use an average length
        #self.update_fish_larvae()

        # Simply move particles with ambient current
        self.advect_ocean_current()

        # Advect particles due to surface wind drag,
        # according to element property wind_drift_factor
        #self.advect_wind()

        # Stokes drift
        #self.stokes_drift()

        # Turbulent Mixing
        if self.get_config('drift:vertical_mixing') is True:
            self.update_terminal_velocity()
            self.vertical_mixing()
#        else:  # Buoyancy
#            self.vertical_buoyancy()

        # Vertical advection
        if self.get_config('drift:vertical_advection') is True:
            self.vertical_advection()

        # Vertical migration
        self.larvae_vertical_migration()

        # Deactivate floats after "drift_days" has passed
        self.deactivate_elements(self.elements.age_seconds > drift_seconds, reason='age > {} days'.format(drift_days))



