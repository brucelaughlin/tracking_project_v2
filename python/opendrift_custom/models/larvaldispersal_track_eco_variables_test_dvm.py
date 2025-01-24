# Turn off vertical migration for now; start with just physics

# Test updated version of dvm, the function for which should also work for staying in the boundary layer (just set target depth)


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

# Setup for random velocity kicks
np.random.seed(0)
#velocity_kick_max_const = 0.2

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
            'drift:max_lifespan_days':{
                'type': 'int',
                'default': 180,
                'min': 1,
                'max': 9999,
                'units': 'days',
                'description': 'Maximum drifter lifespan before deactivation',
                'level': CONFIG_LEVEL_BASIC},
            'drift:random_velocity_kick':{
                'type': 'float',
                'default': 0,
                'min': 0,
                'max': 10,
                'units': 'm/s',
                'description': 'Maximum speed value of random velocity kick, default is 0',
                'level': CONFIG_LEVEL_BASIC},
            'drift:velocity_kick_depth_e_folding_scale':{
                'type': 'float',
                'default': 200,
                'min': 0,
                'max': 1000,
                'units': 'm',
                'description': 'Depth e-folding scale for velocity kicks',
                'level': CONFIG_LEVEL_BASIC},
            #'drift:velocity_kick_depth_max':{
            #    'type': 'float',
            #    'default': 800,
            #    'min': 0,
            #    'max': 1000,
            #    'units': 'm',
            #    'description': 'Maximum depth at which random velocity kicks are applied',
            #    'level': CONFIG_LEVEL_BASIC},
        })

        
        ## IBM configuration options
        #self._add_config({
        #    'IBM:fraction_of_timestep_swimming':
        #        {'type': 'float', 'default': 0.15,
        #         'min': 0.0, 'max': 1.0, 'units': 'fraction',
        #         'description': 'Fraction of timestep swimming',
        #         'level': CONFIG_LEVEL_ADVANCED},
        #    })
        
        #self._set_config_default('drift:vertical_advection', True)
        #self._set_config_default('drift:vertical_mixing', True)
        #self._set_config_default('general:coastline_action', 'previous')
        #self._set_config_default('general:use_auto_landmask', False)

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
    def larvae_vertical_migration(self):

        # Note: hardcoding a guessed average length (1mm) of the larvae; is this how we want to do things?
        larvae_length = 1
        # Note: length is in mm:
        # https://opendrift.github.io/_modules/opendrift/models/larvalfish.html#LarvalFish.larvae_vertical_migration


        # Vertical migration of Larvae
        # Swim function from Peck et al. 2006

        #L = self.elements.length[larvae] # I don't have this implemented, so I hack it:
        L = np.ones(len(self.elements))*larvae_length
       

        # Below, with fraction_of_timestep_swimming = 0.15 (see above), this means
        # that max_migration_per_timestep = 1.75cm... is that too small?

        swim_speed = (0.261*(L**(1.552*L**(-0.08))) - 5.289/L) / 1000
        f = self.get_config('IBM:fraction_of_timestep_swimming')
        max_migration_per_timestep = f*swim_speed*self.time_step.total_seconds()

        # Using here UTC hours. Should be changed to local solar time,
        # although a phase shift of some hours should not make much difference
        if self.time.hour < 12:
            direction = -1  # Swimming down when light is increasing
        else:
            direction = 1  # Swimming up when light is decreasing

        #self.elements.z[larvae] = np.minimum(0, self.elements.z[larvae] + direction*max_migration_per_timestep)
        self.elements.z = np.minimum(0, self.elements.z + direction*max_migration_per_timestep)


    def vertical_migration_taper_function(z_array, target_depth, target_depth_fraction):
        vertical_kicks = 0.5 + 0.5*(np.cos(np.pi/((1-target_depth_fraction) * target_depth) * (z_array - target_depth_fraction*target_depth)))
        return vertical_kicks


    def larvae_vertical_migration(self):
        """Move particles vertically towards target depth according to pre-defined function
        """
        if self.get_config('drift:vertical_kick') == 0:
            logger.debug('Not applying vertical kick towards target depth')
            return


        #local_time = ((self.time.hour-8) + 24) % 24

        # Unlike vertical advection, we want kicks even for particles at the surface (right?)
        #in_ocean = np.where(self.elements.z<0)[0]
        #if len(in_ocean) > 0:

        target_depth_fraction = self.get_config('drift:target_depth_fraction')
        """<target_depth_fraction> is the fraction of the target depth (down from the surface) over which
         the kick should be constant, aside from the random component.  After this fraction of the target
         depth is traversed, the kick velocity tapers off to zero at the target depth,
         according to the chosen function (here, cosine)"""
       
        # Load target depth.  For DVM, this depends on time.  For constant-layer drifting, set day and night values to be equal in config file.
        if self.time.hour < 12:
            target_depth = self.get_config('drift:target_depth_night')
        else:
            target_depth = self.get_config('drift:target_depth_day')


        vkick_max = self.get_config('drift:vertical_kick')
        random_vkick_max = self.get_config('drift:vertical_kick_random')


        # Get the indices of elements in the 4 different regimes specified by the vertical kick function
        upper_top_dex = self.elements.z > -1 * (target_depth * target_depth_fraction) 

        upper_mid_d1 = self.elements.z <= -1 * (target_depth * target_depth_fraction) 
        upper_mid_d2 = self.elements.z > -1 * target_depth
        upper_mid_dex = upper_mid_d1 & upper_mid_d2

        lower_mid_d1 = self.elements.z <= -1 * target_depth
        lower_mid_d2 = self.elements.z > -1 * (target_depth * (2 - target_depth_fraction))
        lower_mid_dex = lower_mid_d1 & lower_mid_d2
        
        lower_bot_dex = self.elements.z <= -1 * (target_depth * (2 - target_depth_fraction))


        # Apply the kicks
        self.elements.z[upper_top_dex] = np.minimum(0,
            (self.elements.z[upper_top_dex] - self.elements.moving[upper_top_dex] 
                * (vkick_max + (np.random.uniform(-1,1,len(self.elements.z[upper_top_dex])) * random_vkick_max)) * self.time_step.total_seconds())

        self.elements.z[upper_mid_dex] = np.minimum(0,
            (self.elements.z[upper_mid_dex] - self.elements.moving[upper_mid_dex] 
                * (vertical_migration_taper_function(self.elements.z[upper_mid_dex]) * vkick_max + (np.random.uniform(-1,1,len(self.elements.z[upper_mid_dex])) * random_vkick_max)) * self.time_step.total_seconds())
            
        self.elements.z[lower_mid_dex] = np.minimum(0,
            (self.elements.z[lower_mid_dex] + self.elements.moving[lower_mid_dex] 
                * (vertical_migration_taper_function(self.elements.z[lower_mid_dex]) * vkick_max + (np.random.uniform(-1,1,len(self.elements.z[lower_mid_dex])) * random_vkick_max)) * self.time_step.total_seconds())
            
        self.elements.z[lower_bot_dex] = np.minimum(0,
            (self.elements.z[lower_bot_dex] + self.elements.moving[lower_bot_dex] 
                * (vkick_max + (np.random.uniform(-1,1,len(self.elements.z[lower_bot_dex])) * random_vkick_max)) * self.time_step.total_seconds())
        
        
        #self.elements.z[upper_top_dex] + self.elements.moving * w * self.time_step.total_seconds())
        

        # From Oceandrift vertical_advection()  method
        #self.elements.z = np.minimum(0,
        #    self.elements.z + self.elements.moving * w * self.time_step.total_seconds())
        
        #else:
        #    logger.debug('No vertical kick for elements at surface')


    def velocity_kick(self):
        # Note that Paul showed me that applying functions (ie square/sqrt) to a uniform function yeilds a pdf which is likely no longer uniform.
        # So, his suggestion was to use a random angle
        
        if self.get_config('drift:random_velocity_kick') == 0:
            logger.debug('Not applying random horizontal velocity tidal kick')
            return


        #kick_mask = np.abs(self.environment.sea_floor_depth_below_sea_level) < self.get_config('drift:velocity_kick_depth_max')
        #kick_mask = kick_mask.astype(int)
       
        #logger.info(f'Min seafloor depth: {np.abs(np.min(self.environment.sea_floor_depth_below_sea_level))}')
        #logger.info(f'Max seafloor depth: {np.abs(np.max(self.environment.sea_floor_depth_below_sea_level))}')

        kick_speeds_pre = self.get_config('drift:random_velocity_kick') / np.exp(np.abs(self.environment.sea_floor_depth_below_sea_level) / self.get_config('drift:velocity_kick_depth_e_folding_scale'))

        kick_speeds = np.random.rand(len(self.elements)) * kick_speeds_pre
        #kick_speeds = kick_mask * np.random.rand(len(self.elements)) * kick_speeds_pre
        #kick_speeds = np.random.rand(len(self.elements)) * self.get_config('drift:random_velocity_kick')
        kick_angles = 2 * np.pi * np.random.rand(len(self.elements))

        x_vel_kicks = np.cos(kick_angles) * kick_speeds
        y_vel_kicks = np.sin(kick_angles) * kick_speeds

        self.update_positions(x_vel_kicks, y_vel_kicks)


    def update(self):
        """Update positions and properties of elements."""
        # copied from "OceanDrift", with the addition of the deactivation

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
        #else:  # Buoyancy
        #    self.vertical_buoyancy()

        # Vertical advection
        if self.get_config('drift:vertical_advection') is True:
            self.vertical_advection()

        # Prescribe velocity kicks
        if self.get_config('drift:random_velocity_kick') > 0:
            self.velocity_kick()

        # How can I print???  Need to confirm that my config settings are being used
        #logger.debug(f"life: {self.get_config('drift:max_lifespan_days')}")
        #logger.debug(f"vadvect: {self.get_config('drift:vertical_advection')}")
        #logger.debug(f"vmix: {self.get_config('drift:vertical_mixing')}")
        #logger.debug(f"coastline action: {self.get_config('general:coastline_action')}")
        #logger.debug(f"internal mask: {self.get_config('general:use_auto_landmask')}")


        # Deactivate floats after "drift_days" has passed
        self.deactivate_elements(self.elements.age_seconds > self.get_config('drift:max_lifespan_days') * 24 * 60 * 60, reason='age > {} days'.format(self.get_config('drift:max_lifespan_days')))
        #self.deactivate_elements(self.elements.age_seconds > drift_seconds, reason='age > {} days'.format(drift_days))


