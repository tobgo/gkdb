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
import json
import datetime
import pandas as pd

db = PostgresqlExtDatabase(database='gkdb', host='gkdb.qualikiz.com')
class BaseModel(Model):
    """A base model that will use our Postgresql database"""
    class Meta:
        database = db
        schema = 'develop'

class Point(BaseModel):
    creator = TextField(help_text='Name of the creator of this entry')
    date = DateTimeField(help_text='Creation date of this entry')
    comment = TextField(help_text='Any comment describing this entry')
    beta = FloatField(help_text='Plasma beta')
    collisionality = FloatField(help_text='Plasma collision frequency')

    def to_dict(self):
        model_dict = model_to_dict(self, exclude=[Point.id, Point.date])
        model_dict['species'] = []
        for species in self.species:
            model_dict['species'].append(
                model_to_dict(species,
                              recurse=False,
                              exclude=[Species.id, Species.point_id]))
        model_dict['wavevectors'] = []
        model_dict['flux_surface'] = model_to_dict(self.flux_surface.get(),
                                                   recurse=False,
                                                   exclude=[Flux_Surface.elongation,
                                                            Flux_Surface.triangularity,
                                                            Flux_Surface.squareness])
        for wavevector in self.wavevector.select():
            model_dict['wavevectors'].append(
                model_to_dict(wavevector,
                              recurse=False,
                              exclude=[Wavevector.id, Wavevector.point_id]))
            eigenvalue_list = model_dict['wavevectors'][-1]['eigenvalues'] = []
            for eigenvalue in wavevector.eigenvalue.select():
                eigenvalue_list.append(
                    model_to_dict(eigenvalue,
                                  recurse=False,
                                  exclude=[Eigenvalue.id,
                                           Eigenvalue.wavevector_id,
                                           Eigenvalue.a_amplitude,
                                           Eigenvalue.phi_amplitude,
                                           Eigenvalue.phi_parity,
                                           Eigenvalue.a_amplitude,
                                           Eigenvalue.a_parity,
                                           Eigenvalue.b_amplitude,
                                           Eigenvalue.b_parity,
                                           ]))
                eigenvalue_list[-1]['eigenvector'] = (
                    model_to_dict(eigenvalue.eigenvector.get(),
                                  recurse=False,
                                  exclude=[Eigenvector.eigenvalue_id]))



        for flux_table in [Particle_Fluxes, Heat_Fluxes_Lab, Heat_Fluxes_Rotating,
                           Momentum_Fluxes_Lab, Momentum_Fluxes_Rotating]:
            name = flux_table.__name__.lower()
            sel = (self.select(Wavevector.id, flux_table)
                       .where(Point.id == self.id)
                       .join(Wavevector, JOIN_LEFT_OUTER)
                       .join(Eigenvalue, JOIN_LEFT_OUTER)
                       .join(Species, JOIN_LEFT_OUTER, (Species.point_id == Point.id))
                       .join(flux_table).tuples())
            if sel.count() > 0:
                model_dict[name] = {}
                df = pd.DataFrame.from_records(list(sel),
                                               columns=['wavevector_id', 'species_id', 'eigenvalue_id',
                                                        'phi_potential', 'a_parallel', 'b_field_parallel'],
                                               index=['wavevector_id', 'species_id', 'eigenvalue_id'])
                xr = df.to_xarray()
                for k, v in xr.data_vars.items():
                    model_dict[name][k] =  v.data.tolist()
            else:
                model_dict[name] = None
        return model_dict

    def to_json(self, path):
        with open(path, 'w') as file_:
            json.dump(self.to_dict(), file_, indent=4, sort_keys=True)

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
    collision_enhancement_factor = FloatField(help_text='Enhancement factor for the collisions of electrons on main ions (to mimic the impact of impurity ions not present in the run)')

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
    toroidal_velocity = FloatField(help_text='Toroidal velocity (common to all species)')
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

def purge_tables():
    clsmembers = inspect.getmembers(sys.modules[__name__], lambda member: inspect.isclass(member) and member.__module__ == __name__)
    for name, cls in clsmembers:
        if name != BaseModel:
            try:
                db.drop_table(cls, cascade=True)
            except ProgrammingError:
                db.rollback()
    db.create_tables([Point, Code, Flux_Surface, Wavevector, Eigenvalue, Eigenvector, Species, Heat_Fluxes_Lab, Momentum_Fluxes_Lab, Heat_Fluxes_Rotating, Momentum_Fluxes_Rotating, Particle_Fluxes])
