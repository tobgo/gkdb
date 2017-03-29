from peewee import *
from peewee import FloatField, FloatField, ProgrammingError
import numpy as np
import inspect
import sys
from playhouse.postgres_ext import PostgresqlExtDatabase, ArrayField
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
    creator = TextField(null=True)
    date = DateTimeField(null=True)
    reference_case_flag = BooleanField(null=True)
    comment = TextField()

class Code(BaseModel):
    point = ForeignKeyField(Point, related_name='code')
    name = TextField(null=True)
    version = TextField(null=True)
    parameters = TextField(null=True)

class Equilibrium(BaseModel):
    point = ForeignKeyField(Point, related_name='equilibrium')
    r_minor = FloatField(null=True)
    q = FloatField(null=True)
    magnetic_shear = FloatField(null=True)
    beta_gradient = FloatField(null=True)
    ip_sign = FloatField(null=True)
    b_field_tor_sign = FloatField(null=True)

class Shape(BaseModel):
    equilibrium = ForeignKeyField(Equilibrium, related_name='shape')
    c0 = FloatField(null=True)
    dc0_dr_minor = FloatField(null=True)
    c = ArrayField(FloatField, null=True)
    s = ArrayField(FloatField, null=True)
    dc_dr_minor = ArrayField(FloatField, null=True)
    ds_dr_minor = ArrayField(FloatField, null=True)

class EmEffects(BaseModel):
    point = ForeignKeyField(Point, related_name='em_effects')
    beta = FloatField(null=True)
    a_parallel_flag = BooleanField(null=True)
    a_perpendicular_flag = BooleanField(null=True)

class Collisions(BaseModel):
    point = ForeignKeyField(Point, related_name='collisions')
    collisionality = FloatField(null=True)
    collisions_enhancement = FloatField(null=True)
    pitch_only_flag = BooleanField(null=True)
    ei_collisions_only_flag = BooleanField(null=True)
    momentum_conservation_flag = BooleanField(null=True)
    energy_conservation_flag = BooleanField(null=True)

class Mode(BaseModel):
    point = ForeignKeyField(Point, related_name='mode')
    radial_wavevector = FloatField(null=True)
    binormal_wavevector = FloatField(null=True)
    number_poloidal_turns = IntegerField(null=True)

class Eigenmode(BaseModel):
    mode                         = ForeignKeyField(Mode, related_name='eigen_mode')
    growth_rate                  = FloatField(null=True)
    frequency                    = FloatField(null=True)
    r_phi_potential_perturbed    = FloatField(null=True) 
    i_phi_potential_perturbed    = FloatField(null=True)
    r_a_parallel_perturbed       = FloatField(null=True)
    i_a_parallel_perturbed       = FloatField(null=True)
    r_b_field_parallel_perturbed = FloatField(null=True)
    i_b_field_parallel_perturbed = FloatField(null=True)
    tolerance_growth_rate        = FloatField(null=True)
    poloidal_angle               = FloatField(null=True)

class Species(BaseModel):
    point = ForeignKeyField(Point, related_name='species')
    charge = FloatField(null=True)
    mass = FloatField(null=True)
    density = FloatField(null=True)
    density_log_gradient = FloatField(null=True)
    temperature = FloatField(null=True)
    temperature_log_gradient = FloatField(null=True)
    toroidal_velocity = FloatField(null=True)
    toroidal_velocity_gradient = FloatField(null=True)

class Fluxes(BaseModel):
    species = ForeignKeyField(Species, related_name='fluxes')
    particle_phi_potential = FloatField(null=True)
    particle_a_parallel = FloatField(null=True)
    particle_b_field_parallel = FloatField(null=True)
    heat_phi_potential = FloatField(null=True)
    heat_a_parallel = FloatField(null=True)
    heat_b_field_parallel = FloatField(null=True)
    kinetic_momentum_phi_potential = FloatField(null=True)
    kinetic_momentum_a_parallel = FloatField(null=True)
    kinetic_momentum_b_field_parallel = FloatField(null=True)
    field_momentum_phi_potential = FloatField(null=True)
    field_momentum_a_parallel = FloatField(null=True)
    field_momentum_b_field_parallel = FloatField(null=True)

def purge_tables():
    clsmembers = inspect.getmembers(sys.modules[__name__], lambda member: inspect.isclass(member) and member.__module__ == __name__)
    for name, cls in clsmembers:
        if name != BaseModel:
            try:
                db.drop_table(cls, cascade=True)
            except ProgrammingError:
                db.rollback()
    db.create_tables([Point, Code, Equilibrium, Shape, EmEffects, Collisions, Mode, Eigenmode, Species, Fluxes])
