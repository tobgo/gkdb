from peewee import *
from peewee import FloatField, FloatField, ProgrammingError
import numpy as np
import inspect
import sys
from playhouse.postgres_ext import PostgresqlExtDatabase, ArrayField, BinaryJSONField
from IPython import embed
import scipy as sc
from scipy import io

db = PostgresqlExtDatabase(database='gkdb', host='gkdb.qualikiz.com')
class BaseModel(Model):
    """A base model that will use our Postgresql database"""
    class Meta:
        database = db
        schema = 'develop'

class Point(BaseModel):
    creator = TextField()
    date = DateTimeField()
    reference_case_flag = BooleanField()
    comment = TextField()
    beta = FloatField()
    collisionality = FloatField()

class Code(BaseModel):
    point = ForeignKeyField(Point, related_name='code')
    name = TextField()
    version = TextField()
    parameters = BinaryJSONField()

    include_centrifugal_effects = BooleanField()
    include_a_parallel = BooleanField()
    include_b_field_parallel = BooleanField()

    collision_pitch_only = BooleanField()
    collision_ei_only = BooleanField()
    collision_momentum_conservation = BooleanField()
    collision_energy_conservation = BooleanField()
    collision_finite_larmor_radius = BooleanField()
    collision_enhancement_factor = FloatField()

    initial_value_run = BooleanField()
    class Meta:
        primary_key = CompositeKey('point')

class Flux_Surface(BaseModel):
    point = ForeignKeyField(Point, related_name='flux_surface')
    r_minor = FloatField()
    # Derived from Shape
    elongation = FloatField()
    triangularity = FloatField()
    squareness = FloatField()
    # Non-derived
    q = FloatField()
    magnetic_shear = FloatField()
    beta_gradient = FloatField()
    ip_sign = SmallIntegerField()
    b_field_tor_sign = SmallIntegerField()
    # Original shape
    c = ArrayField(FloatField)
    s = ArrayField(FloatField)
    dc_dr_minor = ArrayField(FloatField)
    ds_dr_minor = ArrayField(FloatField)
    class Meta:
        primary_key = CompositeKey('point')

class Wavevector(BaseModel):
    point = ForeignKeyField(Point, related_name='wavevector')
    radial_wavevector = FloatField()
    binormal_wavevector = FloatField()
    poloidal_turns = IntegerField()

class Eigenvalue(BaseModel):
    wavevector                   = ForeignKeyField(Wavevector, related_name='eigenvalue')
    growth_rate                  = FloatField()
    frequency                    = FloatField()
    growth_rate_tolerance        = FloatField()

    # Derived quantities
    phi_amplitude = FloatField()
    phi_parity = FloatField()
    a_amplitude = FloatField()
    a_parity = FloatField()
    b_amplitude = FloatField()
    b_parity = FloatField()

class Eigenvector(BaseModel):
    eigenvalue                   = ForeignKeyField(Eigenvalue, related_name='eigenvector')
    r_phi_potential_perturbed    = ArrayField(FloatField)
    i_phi_potential_perturbed    = ArrayField(FloatField)
    r_a_parallel_perturbed       = ArrayField(FloatField, null=True)
    i_a_parallel_perturbed       = ArrayField(FloatField, null=True)
    r_b_field_parallel_perturbed = ArrayField(FloatField, null=True)
    i_b_field_parallel_perturbed = ArrayField(FloatField, null=True)
    poloidal_angle               = ArrayField(FloatField)
    class Meta:
        primary_key = CompositeKey('eigenvalue')

class Species(BaseModel):
    point = ForeignKeyField(Point, related_name='species')
    charge = FloatField()
    mass = FloatField()
    density = FloatField()
    temperature = FloatField()
    toroidal_velocity = FloatField()
    density_log_gradient = FloatField()
    temperature_log_gradient = FloatField()
    toroidal_velocity_gradient = FloatField()

class Particle_Fluxes(BaseModel):
    species = ForeignKeyField(Species, related_name='particle_fluxes')
    eigenvalue = ForeignKeyField(Eigenvalue, related_name='particle_fluxes')
    phi_potential = FloatField()
    a_parallel = FloatField(null=True)
    b_field_parallel = FloatField(null=True)
    class Meta:
        primary_key = CompositeKey('species', 'eigenvalue')

class Heat_Fluxes_Lab(BaseModel):
    species = ForeignKeyField(Species, related_name='heat_fluxes_lab')
    eigenvalue = ForeignKeyField(Eigenvalue, related_name='heat_fluxes_lab')
    phi_potential = FloatField()
    a_parallel = FloatField(null=True)
    b_field_parallel = FloatField(null=True)
    class Meta:
        primary_key = CompositeKey('species', 'eigenvalue')

class Momentum_Fluxes_Lab(BaseModel):
    species = ForeignKeyField(Species, related_name='momentum_fluxes_lab')
    eigenvalue = ForeignKeyField(Eigenvalue, related_name='momentum_fluxes_lab')
    phi_potential = FloatField()
    a_parallel = FloatField(null=True)
    b_field_parallel = FloatField(null=True)
    class Meta:
        primary_key = CompositeKey('species', 'eigenvalue')

class Heat_Fluxes_Rotating(BaseModel):
    species = ForeignKeyField(Species, related_name='heat_fluxes_rotating')
    eigenvalue = ForeignKeyField(Eigenvalue, related_name='heat_fluxes_rotating')
    phi_potential = FloatField()
    a_parallel = FloatField(null=True)
    b_field_parallel = FloatField(null=True)
    class Meta:
        primary_key = CompositeKey('species', 'eigenvalue')

class Momentum_Fluxes_Rotating(BaseModel):
    species = ForeignKeyField(Species, related_name='momentum_fluxes')
    eigenvalue = ForeignKeyField(Eigenvalue, related_name='momentum_fluxes')
    phi_potential = FloatField()
    a_parallel = FloatField(null=True)
    b_field_parallel = FloatField(null=True)
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
purge_tables()
