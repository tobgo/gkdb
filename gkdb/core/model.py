from peewee import *
from peewee import DoubleField, FloatField, ProgrammingError
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
    output_flag = BooleanField(null=True)

class Equilibrium(BaseModel):
    point = ForeignKeyField(Point, related_name='equilibrium')
    r_minor = DoubleField(null=True)
    q = DoubleField(null=True)
    magnetic_shear = DoubleField(null=True)
    beta_gradient = DoubleField(null=True)
    ip_sign = FloatField(null=True)
    b_field_tor_sign = FloatField(null=True)

class Shape(BaseModel):
    equilibrium = ForeignKeyField(Equilibrium, related_name='shape')
    c0 = DoubleField(null=True)
    dc0_dr_minor = DoubleField(null=True)
    c = ArrayField(DoubleField, null=True)
    s = ArrayField(DoubleField, null=True)
    dc_dr_minor = ArrayField(DoubleField, null=True)
    ds_dr_minor = ArrayField(DoubleField, null=True)

class EmEffects(BaseModel):
    point = ForeignKeyField(Point, related_name='em_effects')
    beta = DoubleField(null=True)
    a_parallel_flag = BooleanField(null=True)
    a_perpendicular_flag = BooleanField(null=True)

class Collisions(BaseModel):
    point = ForeignKeyField(Point, related_name='collisions')
    collisionality = DoubleField(null=True)
    collisions_enhancement = DoubleField(null=True)
    pitch_only_flag = BooleanField(null=True)
    ei_collisions_only_flag = BooleanField(null=True)
    momentum_conservation_flag = BooleanField(null=True)
    energy_conservation_flag = BooleanField(null=True)

class MasterMode(BaseModel):
    point = ForeignKeyField(Point, related_name='master_mode')
    radial_wavevector = DoubleField(null=True)
    binormal_wavevector = DoubleField(null=True)

class Mode(BaseModel):
    mastermode                   = ForeignKeyField(MasterMode, related_name='modes')
    growth_rate                  = DoubleField(null=True)
    frequency                    = DoubleField(null=True)
    r_phi_potential_perturbed    = DoubleField(null=True) 
    i_phi_potential_perturbed    = DoubleField(null=True)
    r_a_parallel_perturbed       = DoubleField(null=True)
    i_a_parallel_perturbed       = DoubleField(null=True)
    r_b_field_parallel_perturbed = DoubleField(null=True)
    i_b_field_parallel_perturbed = DoubleField(null=True)

class Species(BaseModel):
    point = ForeignKeyField(Point, related_name='species')
    # Input
    charge = DoubleField(null=True)
    mass = DoubleField(null=True)
    density = DoubleField(null=True)
    density_log_gradient = DoubleField(null=True)
    temperature = DoubleField(null=True)
    temperature_log_gradient = DoubleField(null=True)
    toroidal_velocity = DoubleField(null=True)
    toroidal_velocity_gradient = DoubleField(null=True)
    # Output
    particle_phi_potential = DoubleField(null=True)
    particle_a_parallel = DoubleField(null=True)
    particle_b_field_parallel = DoubleField(null=True)
    heat_phi_potential = DoubleField(null=True)
    heat_a_parallel = DoubleField(null=True)
    heat_b_field_parallel = DoubleField(null=True)
    kinetic_momentum_phi_potential = DoubleField(null=True)
    kinetic_momentum_a_parallel = DoubleField(null=True)
    kinetic_momentum_b_field_parallel = DoubleField(null=True)
    field_momentum_phi_potential = DoubleField(null=True)
    field_momentum_a_parallel = DoubleField(null=True)
    field_momentum_b_field_parallel = DoubleField(null=True)

def purge_tables():
    clsmembers = inspect.getmembers(sys.modules[__name__], lambda member: inspect.isclass(member) and member.__module__ == __name__)
    for name, cls in clsmembers:
        if name != BaseModel:
            try:
                db.drop_table(cls, cascade=True)
            except ProgrammingError:
                db.rollback()
    db.create_tables([Point, Code, Equilibrium, Shape, EmEffects, Collisions, MasterMode, Mode, Species])
