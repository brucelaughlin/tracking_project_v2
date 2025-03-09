# Turn off vertical migration for now; start with just physics

# Test updated version of dvm, the function for which should also work for staying in the boundary layer (just set target depth)


#import datetime
#import numpy as np
#import logging; logger = logging.getLogger(__name__)
#from opendrift.models.oceandrift import Lagrangian3DArray, OceanDrift
#from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED

import sys
import datetime
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
np.random.seed(int(datetime.datetime.now().timestamp()))

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
            'drift:horizontal_velocity_kick_max':{
                'type': 'float',
                'default': 0,
                'min': 0,
                'max': 1,
                'units': 'm/s',
                'description': 'Maximum speed value of random velocity kick, default is 0',
                'level': CONFIG_LEVEL_BASIC},
            'drift:horizontal_velocity_kick_depth_e_folding_scale':{
                'type': 'float',
                'default': 200,
                'min': 0,
                'max': 1000,
                'units': 'm',
                'description': 'Depth e-folding scale for velocity kicks',
                'level': CONFIG_LEVEL_BASIC},
            'drift:vertical_swim_speed':{
                'type': 'float',
                'default': 0,
                'min': 0,
                'max': 1,
                'units': 'm/s',
                'description': 'Vertical velocity kick speed, default is 0',
                'level': CONFIG_LEVEL_BASIC},
            'drift:vertical_migration_random':{
                'type': 'float',
                'default': 0,
                'min': 0,
                'max': 1,
                'units': 'm/s',
                'description': 'Random error to add to vertical velocity kick speed, default is 0',
                'level': CONFIG_LEVEL_BASIC},
            'drift:target_depth_day':{
                'type': 'float',
                'default': 0,
                'min': 0,
                'max': 1000,
                'units': 'm/s',
                'description': 'Daytime target depth for larvae, default is 0',
                'level': CONFIG_LEVEL_BASIC},
            'drift:target_depth_night':{
                'type': 'float',
                'default': 0,
                'min': 0,
                'max': 1000,
                'units': 'm/s',
                'description': 'Nighttime target depth for larvae, default is 0',
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

        
    
    def vertical_migration_function(self, z_array, target_depth):
        swim_cut_distance_from_target_depth = 5
        kick_cutoff_h_minus_target_depth = -5

        vertical_swim_speed_max = self.get_config('drift:vertical_swim_speed')
        timestep_seconds = self.time_step.total_seconds()
        sign_mask = np.ones(len(z_array))
        sign_mask[z_array < target_depth] *= -1

        upper_swim_cut_depth = target_depth - swim_cut_distance_from_target_depth
        lower_swim_cut_depth = target_depth + swim_cut_distance_from_target_depth
        
        random_vertical_kick_std_timestep = swim_cut_distance_from_target_depth/2

        swim_tanh = vertical_swim_speed_max * (-1 * (1 - (1 + np.tanh(z_array - upper_swim_cut_depth))/2  - (1 + np.tanh(z_array - lower_swim_cut_depth))/2))
        swim_line = (z_array - target_depth)/timestep_seconds
        swim_velocity_use = np.minimum(abs(swim_tanh), abs(swim_line)) * sign_mask

        swim_velocity_random_kick = np.random.normal(0, random_vertical_kick_std_timestep/timestep_seconds, len(z_array))
        
        # Prevent vertical kicks if target depth is more than abs(kick_cutoff_h_minus_target_depth) meters below seafloor 
        no_vertical_kick_mask = (self.environment.sea_floor_depth_below_sea_level - target_depth) < kick_cutoff_h_minus_target_depth
        swim_velocity_random_kick = swim_velocity_random_kick * np.logical_not(no_vertical_kick_mask)

        #vertical_displacements = swim_velocity_use * timestep_seconds
        vertical_displacements = (swim_velocity_use + swim_velocity_random_kick) * timestep_seconds

        
        vertical_displacement_dict = dict()
        vertical_displacement_dict['swim_velocities'] = swim_velocity_use
        vertical_displacement_dict['random_kick'] = swim_velocity_random_kick


        return vertical_displacement_dict
        #return vertical_displacements
        

    def vertical_swimming(self):
        """Move particles vertically towards target depth according to pre-defined function
        """
        #if self.get_config('drift:vertical_migration') == 0:
        #    logger.debug('Not applying vertical kick towards target depth')
        #    return

        # Load target depth.  For DVM, this depends on time.  For constant-layer drifting, set day and night values to be equal in config file.
        if self.time.hour < 6 or self.time.hour >= 18:
            target_depth = self.get_config('drift:target_depth_night')
        else:
            target_depth = self.get_config('drift:target_depth_day')

        #random_vertical_swim_speed_max = self.get_config('drift:vertical_migration_random')

        # The vertical kick function works with positive depths, so multiply the negative depths used by Opendrifty by -1 during passing of argument
        #vertical_displacements = self.vertical_migration_function(-1 * self.elements.z, target_depth)
        
        vertical_displacement_dict = self.vertical_migration_function(-1 * self.elements.z, target_depth)
       
        # Update vertical position in two steps.  First, apply swimming.  Then, the random kick.  In this way, we allow kicks even if particles have gotten stuck
        # on the bottom by the Opendrift (move to seafloor) algorithm.  

        self.elements.z = np.minimum(0, self.elements.z + vertical_displacement_dict['swim_velocities'] * self.time_step.total_seconds())
       
        # Create boolean mask to flip negative kicks if a particle is stuck on the bottom (my thinking - there's no point in kicking downwards if they're already at the bottm.)
        # Maybe this isn't too statistically defensible, but it seems better to have a cloud above the floor than a half cloud, half stuck particles

        stuck_on_bottom_mask = self.elements.z == self.environment.sea_floor_depth_below_sea_level
        vertical_displacement_dict['random_kick'][stuck_on_bottom_mask] = abs(vertical_displacement_dict['random_kick'][stuck_on_bottom_mask])
        
        self.elements.z = np.minimum(0, self.elements.z + vertical_displacement_dict['random_kick'] * self.time_step.total_seconds())
        
        logger.info(f'Seafloor depths: {self.environment.sea_floor_depth_below_sea_level}')
        
        #self.elements.z = np.minimum(0, self.elements.z + vertical_displacements)

    
        
    def horizontal_velocity_kick(self):
        # Note that Paul showed me that applying functions (ie square/sqrt) to a uniform function yeilds a pdf which is likely no longer uniform.
        # So, his suggestion was to use a random angle
        
        #if self.get_config('drift:horizontal_velocity_kick_max') == 0:
        #    logger.debug('Not applying random horizontal velocity tidal kick')
        #    return


        #kick_mask = np.abs(self.environment.sea_floor_depth_below_sea_level) < self.get_config('drift:velocity_kick_depth_max')
        #kick_mask = kick_mask.astype(int)
       
        #logger.info(f'Min seafloor depth: {np.abs(np.min(self.environment.sea_floor_depth_below_sea_level))}')
        #logger.info(f'Max seafloor depth: {np.abs(np.max(self.environment.sea_floor_depth_below_sea_level))}')

        kick_speeds_pre = self.get_config('drift:horizontal_velocity_kick_max') / np.exp(np.abs(self.environment.sea_floor_depth_below_sea_level) / self.get_config('drift:horizontal_velocity_kick_depth_e_folding_scale'))

        kick_speeds = np.random.rand(len(self.elements)) * kick_speeds_pre
        #kick_speeds = kick_mask * np.random.rand(len(self.elements)) * kick_speeds_pre
        #kick_speeds = np.random.rand(len(self.elements)) * self.get_config('drift:horizontal_velocity_kick_max')
        kick_angles = 2 * np.pi * np.random.rand(len(self.elements))

        x_vel_kicks = np.cos(kick_angles) * kick_speeds
        y_vel_kicks = np.sin(kick_angles) * kick_speeds

        self.update_positions(x_vel_kicks, y_vel_kicks)


    def update(self):
        """Update positions and properties of elements."""
        # copied from "OceanDrift", with the addition of the deactivation

        # Simply move particles with ambient current
        self.advect_ocean_current()

        # Stokes drift
        #self.stokes_drift()


        # Vertical advection
        if self.get_config('drift:vertical_advection') is True:
            self.vertical_advection()

        # Prescribe horizontal velocity kicks
        if self.get_config('drift:horizontal_velocity_kick_max') > 0:
            self.horizontal_velocity_kick()
        
        # Prescribe vertical velocity kicks
        if self.get_config('drift:vertical_swim_speed') > 0:
            self.vertical_swimming()
        

        # Turbulent Mixing
        if self.get_config('drift:vertical_mixing') is True:
            self.update_terminal_velocity()
            self.vertical_mixing()
        #else:  # Buoyancy
        #    self.vertical_buoyancy()

        # Attempting to print memory usage to log files
        memory_usage_current = self.memory_usage[-1]
        logger.info(f"       MEMORY USAGE (GB)      : {memory_usage_current}")


        # How can I print???  Need to confirm that my config settings are being used
        #logger.debug(f"life: {self.get_config('drift:max_lifespan_days')}")
        #logger.debug(f"vadvect: {self.get_config('drift:vertical_advection')}")
        #logger.debug(f"vmix: {self.get_config('drift:vertical_mixing')}")
        #logger.debug(f"coastline action: {self.get_config('general:coastline_action')}")
        #logger.debug(f"internal mask: {self.get_config('general:use_auto_landmask')}")


        # Deactivate floats after "drift_days" has passed
        self.deactivate_elements(self.elements.age_seconds > self.get_config('drift:max_lifespan_days') * 24 * 60 * 60, reason='age > {} days'.format(self.get_config('drift:max_lifespan_days')))
        #self.deactivate_elements(self.elements.age_seconds > drift_seconds, reason='age > {} days'.format(drift_days))


