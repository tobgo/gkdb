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

class Code(BaseModel):
    name = TextField()
    version = TextField()
    parameters = BinaryJSONField()

class Flux_Surface(BaseModel):
    r_minor = FloatField()
    q = FloatField()
    magnetic_shear = FloatField()
    beta_gradient = FloatField()
    ip_sign = SmallIntegerField()
    b_field_tor_sign = SmallIntegerField()

class Shape(BaseModel):
    equilibrium = ForeignKeyField(Flux_Surface, related_name='shape')
    c = ArrayField(FloatField)
    s = ArrayField(FloatField)
    dc_dr_minor = ArrayField(FloatField)
    ds_dr_minor = ArrayField(FloatField)

class Em_Effects(BaseModel):
    beta = FloatField()
    a_parallel_flag = BooleanField()
    b_field_parallel_flag = BooleanField()

class Collisions(BaseModel):
    collisionality = FloatField()
    collisions_enhancement = FloatField()
    pitch_only_flag = BooleanField()
    ei_collisions_only_flag = BooleanField()
    momentum_conservation_flag = BooleanField()
    energy_conservation_flag = BooleanField()

class Mode(BaseModel):
    radial_wavevector = FloatField()
    binormal_wavevector = FloatField()
    number_poloidal_turns = IntegerField()

class Eigenmode(BaseModel):
    mode                         = ForeignKeyField(Mode, related_name='eigen_mode')
    growth_rate                  = FloatField()
    frequency                    = FloatField()
    tolerance_growth_rate        = FloatField()
    r_phi_potential_perturbed    = ArrayField(FloatField)
    i_phi_potential_perturbed    = ArrayField(FloatField)
    r_a_parallel_perturbed       = ArrayField(FloatField, null=True)
    i_a_parallel_perturbed       = ArrayField(FloatField, null=True)
    r_b_field_parallel_perturbed = ArrayField(FloatField, null=True)
    i_b_field_parallel_perturbed = ArrayField(FloatField, null=True)
    poloidal_angle               = ArrayField(FloatField)

class Species(BaseModel):
    charge = FloatField()
    mass = FloatField()
    density = FloatField()
    temperature = FloatField()
    toroidal_velocity = FloatField()

class Gradients(BaseModel):
    density_log_gradient = FloatField()
    temperature_log_gradient = FloatField()
    toroidal_velocity_gradient = FloatField()

class Fluxes(BaseModel):
    species = ForeignKeyField(Species, related_name='fluxes')
    mode = ForeignKeyField(Mode, related_name='fluxes')
    particle_phi_potential = FloatField()
    particle_a_parallel = FloatField(null=True)
    particle_b_field_parallel = FloatField(null=True)
    heat_phi_potential = FloatField()
    heat_a_parallel = FloatField(null=True)
    heat_b_field_parallel = FloatField(null=True)
    kinetic_momentum_phi_potential = FloatField()
    kinetic_momentum_a_parallel = FloatField(null=True)
    kinetic_momentum_b_field_parallel = FloatField(null=True)

class Point(BaseModel):
    creator = TextField()
    date = DateTimeField()
    reference_case_flag = BooleanField()
    comment = TextField()
    related_name = 'point'
    code = ForeignKeyField(Code, related_name=related_name)
    flux_surface = ForeignKeyField(Flux_Surface, related_name=related_name)
    em_effects = ForeignKeyField(Em_Effects, related_name=related_name)
    collisions = ForeignKeyField(Collisions, related_name=related_name)
    mode = ForeignKeyField(Mode, related_name=related_name)
    species = ForeignKeyField(Species, related_name=related_name)
    gradients = ForeignKeyField(Gradients, related_name=related_name)


def purge_tables():
    clsmembers = inspect.getmembers(sys.modules[__name__], lambda member: inspect.isclass(member) and member.__module__ == __name__)
    for name, cls in clsmembers:
        if name != BaseModel:
            try:
                db.drop_table(cls, cascade=True)
            except ProgrammingError:
                db.rollback()
    db.create_tables([Point, Code, Flux_Surface, Shape, Em_Effects, Collisions, Mode, Eigenmode, Species, Fluxes, Gradients])
purge_tables()
