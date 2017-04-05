from peewee import *
from peewee import FloatField, FloatField, ProgrammingError
import numpy as np
import inspect
import sys
from playhouse.postgres_ext import PostgresqlExtDatabase, ArrayField, BinaryJSONField
from playhouse.shortcuts import model_to_dict
from IPython import embed
import scipy as sc
from scipy import io

db = PostgresqlExtDatabase(database='gkdb', host='gkdb.qualikiz.com')
class BaseModel(Model):
    """A base model that will use our Postgresql database"""
    class Meta:
        database = db
        schema = 'develop'

class Tag(BaseModel):
    name =    TextField(null=True)
    doi =     TextField(null=True)
    comment = TextField(null=True)

class Point(BaseModel):
    creator = TextField(help_text='Name of the creator of this entry')
    date = DateTimeField(help_text='Creation date of this entry')
    comment = TextField(help_text='Any comment describing this entry')

class Point_Tag(BaseModel):
    point = ForeignKeyField(Point)
    tag = ForeignKeyField(Tag)

class Code(BaseModel):
    point = ForeignKeyField(Point, related_name='code')
    name = TextField(help_text='Name of the code used for this entry')
    version = TextField(help_text='Version of the code used for this entry')
    parameters = BinaryJSONField(help_text='Key/value store containing the code dependent inputs')

    include_centrifugal_effects = BooleanField(help_text='True if centrifugal effects were included, false otherwise')
    include_a_parallel = BooleanField(help_text='True if fluctuations of the parallel vector potential (magnetic flutter) were retained, false otherwise')
    include_b_field_parallel = BooleanField(help_text='True if fluctuations of the parallel magnetic field (magnetic compression) were retained, false otherwise')

    collision_pitch_only = BooleanField(help_text='True if pitch-angle scattering only was retained in the collision operator, false otherwise')
    collision_ei_only = BooleanField(help_text='True if electron to main ion collisions only were retained in the collision operator. False if all species collisions were retained')
    collision_momentum_conservation = BooleanField(help_text='True if the collision operator conserves momentum, false otherwise.')
    collision_energy_conservation = BooleanField(help_text='True if the collision operator conserves energy, false otherwise.')
    collision_finite_larmor_radius = BooleanField(help_text='True if the collision operator includes finite Larmor radius effects, false otherwise.')

    initial_value_run = BooleanField(help_text='True if the run was an initial value run. False if it was an eigenvalue run.')
    class Meta:
        primary_key = CompositeKey('point')

class Flux_Surface(BaseModel):
    point = ForeignKeyField(Point, related_name='flux_surface')
    r_minor = FloatField(help_text='Minor radius of the flux surface of interest')
    # Derived from Shape
    elongation = FloatField(help_text='Elongation of the flux surface of interest. Computed internally from the shape parameters (c_n,s_n)')
    triangularity = FloatField(help_text='Triangularity of the flux surface of interest. Computed internally from the shape parameters (c_n,s_n)')
    squareness = FloatField(help_text='Squareness of the flux surface of interest. Computed internally from the shape parameters (c_n,s_n)')
    # Non-derived
    q = FloatField(help_text='Safety factor')
    magnetic_shear = FloatField(help_text='Magnetic shear')
    pressure_gradient = FloatField(help_text='Total pressure gradient (with respect to r_minor) used to characterise the local magnetic equilibrium')
    ip_sign = SmallIntegerField(help_text='Direction of the toroidal plasma current, positive when anticlockwise from above')
    b_field_tor_sign = SmallIntegerField(help_text='Direction of the toroidal magnetic field, positive when anticlockwise from above')
    # Original shape
    c = ArrayField(FloatField, help_text='Array containing the c_n coefficients parametrising the flux surface of interest. ')
    s = ArrayField(FloatField, help_text='Array containing the s_n coefficients parametrising the flux surface of interest. The first element is always zero.')
    dc_dr_minor = ArrayField(FloatField, help_text='Radial derivative (with respect to r_minor) of the c_n coefficients')
    ds_dr_minor = ArrayField(FloatField, help_text='Radial derivative (with respect to r_minor) of the s_n coefficients. The first element is always zero.')
    class Meta:
        primary_key = CompositeKey('point')

class Species_Global(BaseModel):
    point = ForeignKeyField(Point, related_name='species_global')
    beta = FloatField(help_text='Plasma beta')
    collisionality = FloatField(help_text='Plasma collision frequency')
    collision_enhancement_factor = FloatField(help_text='Enhancement factor for the collisions of electrons on main ions (to mimic the impact of impurity ions not present in the run)')
    toroidal_velocity = FloatField(help_text='Toroidal velocity (common to all species)')
    # Derived from Species
    zeff = FloatField(null=True)
    class Meta:
        primary_key = CompositeKey('point')

class Wavevector(BaseModel):
    point = ForeignKeyField(Point, related_name='wavevector')
    radial_wavevector = FloatField(help_text='Radial component of the wavevector')
    binormal_wavevector = FloatField(help_text='Binormal component of the wavevector')
    poloidal_turns = IntegerField(help_text='Number of poloidal turns covered by the flux-tube domain (i.e. number of coupled radial modes included in the simulation)')

class Eigenvalue(BaseModel):
    wavevector                   = ForeignKeyField(Wavevector, related_name='eigenvalue')
    growth_rate                  = FloatField(help_text='Mode growth rate')
    frequency                    = FloatField(help_text='Mode frequency')
    growth_rate_tolerance        = FloatField(help_text='Tolerance used for to determine the mode growth rate convergence')

    # Derived quantities
    phi_amplitude = FloatField(help_text='Relative amplitude of the electrostatic potential perturbations. Computed internally from the parallel structure of the fields (eigenvectors)')
    phi_parity = FloatField(help_text='Parity of the electrostatic potential perturbations. Computed internally from the parallel structure of the fields (eigenvectors)')
    a_amplitude = FloatField(help_text='Relative amplitude of the parallel vector potential perturbations. Computed internally from the parallel structure of the fields (eigenvectors)')
    a_parity = FloatField(help_text='Parity of the parallel vector potential perturbations. Computed internally from the parallel structure of the fields (eigenvectors)')
    b_amplitude = FloatField(help_text='Relative amplitude of the parallel magnetic field perturbations. Computed internally from the parallel structure of the fields (eigenvectors)')
    b_parity = FloatField(help_text='Parity of the parallel magnetic field perturbations. Computed internally from the parallel structure of the fields (eigenvectors)')

class Eigenvector(BaseModel):
    eigenvalue                   = ForeignKeyField(Eigenvalue, related_name='eigenvector')
    r_phi_potential_perturbed    = ArrayField(FloatField, help_text='Parallel structure of the electrostatic potential perturbations (real part)')
    i_phi_potential_perturbed    = ArrayField(FloatField, help_text='Parallel structure of the electrostatic potential perturbations (imaginary part)')
    r_a_parallel_perturbed       = ArrayField(FloatField, null=True, help_text='Parallel structure of the parallel vector potential perturbations (real part)')
    i_a_parallel_perturbed       = ArrayField(FloatField, null=True, help_text='Parallel structure of the parallel vector potential perturbations (imaginary part)')
    r_b_field_parallel_perturbed = ArrayField(FloatField, null=True, help_text='Parallel structure of the parallel magnetic field perturbations (real part)')
    i_b_field_parallel_perturbed = ArrayField(FloatField, null=True, help_text='Parallel structure of the parallel magnetic field perturbations (imaginary part)')
    poloidal_angle               = ArrayField(FloatField, help_text='Poloidal angle grid used to specify the parallel structure of the fields (eigenvectors)')
    class Meta:
        primary_key = CompositeKey('eigenvalue')

class Species(BaseModel):
    point = ForeignKeyField(Point, related_name='species')
    charge = FloatField(help_text='Species charge')
    mass = FloatField(help_text='Species mass')
    density = FloatField(help_text='Species density')
    temperature = FloatField(help_text='Species temperature')
    density_log_gradient = FloatField(help_text='Species logarithmic density gradient (with respect to r_minor)')
    temperature_log_gradient = FloatField(help_text='Species logarithmic temperature gradient (with respect to r_minor)')
    toroidal_velocity_gradient = FloatField(help_text='Species toroidal velocity gradient (with respect to r_minor)')

class Particle_Fluxes(BaseModel):
    species = ForeignKeyField(Species, related_name='particle_fluxes')
    eigenvalue = ForeignKeyField(Eigenvalue, related_name='particle_fluxes')
    phi_potential = FloatField(help_text='Gyrocenter particle flux due to the electrostatic potential fluctuations. Identical in the Laboratory and rotating frames')
    a_parallel = FloatField(null=True, help_text='Gyrocenter particle flux due to the parallel vector potential fluctuations (magnetic flutter). Identical in the Laboratory and rotating frames')
    b_field_parallel = FloatField(null=True, help_text='Gyrocenter particle flux due to the parallel magnetic field fluctuations (magnetic compression). Identical in the Laboratory and rotating frames')
    class Meta:
        primary_key = CompositeKey('species', 'eigenvalue')

class Heat_Fluxes_Lab(BaseModel):
    species = ForeignKeyField(Species, related_name='heat_fluxes_lab')
    eigenvalue = ForeignKeyField(Eigenvalue, related_name='heat_fluxes_lab')
    phi_potential = FloatField(help_text='Gyrocenter heat flux due to the electrostatic potential fluctuations in the Laboratory frame')
    a_parallel = FloatField(null=True, help_text='Gyrocenter heat flux due to the parallel vector potential fluctuations (magnetic flutter) in the Laboratory frame')
    b_field_parallel = FloatField(null=True, help_text='Gyrocenter heat flux due to the parallel magnetic field fluctuations (magnetic compression) in the Laboratory frame')
    class Meta:
        primary_key = CompositeKey('species', 'eigenvalue')

class Momentum_Fluxes_Lab(BaseModel):
    species = ForeignKeyField(Species, related_name='momentum_fluxes_lab')
    eigenvalue = ForeignKeyField(Eigenvalue, related_name='momentum_fluxes_lab')
    phi_potential = FloatField(help_text='Gyrocenter momentum flux (toroidal projection of the parallel contribution only) due to the electrostatic potential fluctuations in the Laboratory frame')
    a_parallel = FloatField(null=True, help_text='Gyrocenter momentm flux (toroidal projection of the parallel contribution only) due to the parallel vector potential fluctuations (magnetic flutter) in the Laboratory frame')
    b_field_parallel = FloatField(null=True, help_text='Gyrocenter momentum flux (toroidal projection of the parallel contribution only) due to the parallel magnetic field fluctuations (magnetic compression) in the Laboratory frame')
    class Meta:
        primary_key = CompositeKey('species', 'eigenvalue')

class Heat_Fluxes_Rotating(BaseModel):
    species = ForeignKeyField(Species, related_name='heat_fluxes_rotating')
    eigenvalue = ForeignKeyField(Eigenvalue, related_name='heat_fluxes_rotating')
    phi_potential = FloatField(help_text='Gyrocenter heat flux due to the electrostatic potential fluctuations in the rotating frame')
    a_parallel = FloatField(null=True, help_text='Gyrocenter heat flux due to the parallel vector potential fluctuations (magnetic flutter) in the rotating frame')
    b_field_parallel = FloatField(null=True, help_text='Gyrocenter heat flux due to the parallel magnetic field fluctuations (magnetic compression) in the rotating frame')
    class Meta:
        primary_key = CompositeKey('species', 'eigenvalue')

class Momentum_Fluxes_Rotating(BaseModel):
    species = ForeignKeyField(Species, related_name='momentum_fluxes')
    eigenvalue = ForeignKeyField(Eigenvalue, related_name='momentum_fluxes')
    phi_potential = FloatField(help_text='Gyrocenter momentum flux (toroidal projection of the parallel contribution only) due to the electrostatic potential fluctuations in the rotating frame')
    a_parallel = FloatField(null=True, help_text='Gyrocenter momentm flux (toroidal projection of the parallel contribution only) due to the parallel vector potential fluctuations (magnetic flutter) in the rotating frame')
    b_field_parallel = FloatField(null=True, help_text='Gyrocenter momentum flux (toroidal projection of the parallel contribution only) due to the parallel magnetic field fluctuations (magnetic compression) in the rotating frame')
    class Meta:
        primary_key = CompositeKey('species', 'eigenvalue')

class Moments_Rotating(BaseModel):
    species = ForeignKeyField(Species, related_name='moments')
    eigenvalue = ForeignKeyField(Eigenvalue, related_name='moments')
    r_density                      = ArrayField(FloatField)
    r_parallel_velocity            = ArrayField(FloatField)
    r_parallel_temperature         = ArrayField(FloatField)
    r_perpendicular_temperature    = ArrayField(FloatField)
    r_j0_density                   = ArrayField(FloatField)
    r_j0_parallel_velocity         = ArrayField(FloatField)
    r_j0_parallel_temperature      = ArrayField(FloatField)
    r_j0_perpendicular_temperature = ArrayField(FloatField)
    i_density                      = ArrayField(FloatField)
    i_parallel_velocity            = ArrayField(FloatField)
    i_parallel_temperature         = ArrayField(FloatField)
    i_perpendicular_temperature    = ArrayField(FloatField)
    i_j0_density                   = ArrayField(FloatField)
    i_j0_parallel_velocity         = ArrayField(FloatField)
    i_j0_parallel_temperature      = ArrayField(FloatField)
    i_j0_perpendicular_temperature = ArrayField(FloatField)
    class Meta:
        primary_key = CompositeKey('species', 'eigenvalue')

def purge_tables():
    clsmembers = inspect.getmembers(sys.modules[__name__], lambda member: inspect.isclass(member) and member.__module__ == __name__)
    for name, cls in clsmembers:
        if name != BaseModel:
            try:
                db.drop_table(cls, cascade=True)
            except ProgrammingError:
                db.rollback()
    db.create_tables([Tag, Point_Tag, Point, Code, Flux_Surface, Wavevector, Eigenvalue, Eigenvector, Species, Heat_Fluxes_Lab, Momentum_Fluxes_Lab, Heat_Fluxes_Rotating, Momentum_Fluxes_Rotating, Particle_Fluxes, Moments_Rotating, Species_Global])
