from peewee import *
from peewee import DoubleField, FloatField, ProgrammingError
import numpy as np
import inspect
import sys
from playhouse.postgres_ext import PostgresqlExtDatabase, ArrayField
from IPython import embed
import scipy as sc
from scipy import io
import csv
import os
sys.path.append(os.path.dirname(os.path.realpath(os.path.join(__file__, '../..'))))
from gkdb.core.model import *
def flatten_floatstruct(floatstruct):
    nparray = []
    for ii, el in enumerate(floatstruct):
        if float(el) == -9e40 or float(el) == -999999999:
            val = None
        else:
            val = float(el)
        nparray.append(val)
    return nparray

def matdict_to_SQL(matdict, eigenfunc_line):
    if np.all(matdict[3] == 0):
        caseflag = False
    elif np.all(matdict[3] == 1):
        caseflag = True
    else:
        caseflag = None

    if str(matdict[0][0][0][0]) == '[]':
        comment = None
    else:
        raise NotImplementedError

    try:
        creator = matdict[1][0]
    except IndexError:
        creator = None
    try:
        date = matdict[2][0]
    except IndexError:
        date = None
    point = Point(comment=matdict[0], creator=creator, date=date, reference_case_flag=caseflag)
    point.save()

    inputs = matdict[4][0][0]
    equi = inputs[0][0][0]
    sha = equi[4][0][0]
    equi = flatten_floatstruct([equi[0], equi[1], equi[2], equi[3], equi[5], equi[6]])

    equilibrium = Equilibrium(point=point, r_minor=equi[0], q=equi[1], magnetic_shear=equi[2], beta_gradient=equi[3], ip_sign=equi[4], b_field_tor_sign=equi[5])
    equilibrium.save()
    npsha = flatten_floatstruct(sha)


    shape = Shape(equilibrium=equilibrium, c0=npsha[0], dc0_dr_minor=npsha[1], c=npsha[2::4], s=npsha[3::4], dc_dr_minor=npsha[4::4], ds_dr_minor=npsha[5::4])
    shape.save()

    spe = inputs[1][0][0]
    npspe = flatten_floatstruct(spe)

    em_effects = inputs[2][0][0]
    em_effects = flatten_floatstruct(em_effects)
    em_effects = EmEffects(point=point, beta=em_effects[0], a_parallel_flag=em_effects[1], a_perpendicular_flag=em_effects[2])
    em_effects.save()

    
    collisions = inputs[3][0][0]
    collisions = flatten_floatstruct(collisions)
    collisions = Collisions(point=point,
                            collisionality            =collisions[0],
                            collisions_enhancement    =collisions[1],
                            pitch_only_flag           =collisions[0],
                            ei_collisions_only_flag   =collisions[1],
                            momentum_conservation_flag=collisions[2],
                            energy_conservation_flag  =collisions[3])
    collisions.save()

    mode = inputs[4][0][0]
    mode = flatten_floatstruct(mode)
    mastermode = MasterMode(point=point, radial_wavevector=mode[0], binormal_wavevector=mode[1])
    mastermode.save()
    outputs = matdict[5][0][0]
    eigenvalues = flatten_floatstruct(outputs[0][0][0])
    eigenfunc = flatten_floatstruct(eigenfunc_line[2:])
    if eigenfunc[1] == -9e40:
        poloidal_angle = None
    else:
        poloidal_angle = float(eigenfunc[1])


    for val in np.hstack([np.array(eigenvalues).reshape(3,2), np.array(eigenfunc).reshape(3,6)]):
        mode = Mode(mastermode=mastermode, growth_rate=val[0], frequency=val[1], poloidal_angle=poloidal_angle,
                    r_phi_potential_perturbed   = val[2], 
                    i_phi_potential_perturbed   = val[3],
                    r_a_parallel_perturbed      = val[4],
                    i_a_parallel_perturbed      = val[5],
                    r_b_field_parallel_perturbed= val[6],
                    i_b_field_parallel_perturbed= val[7])
        mode.save()
    eigenvalues = flatten_floatstruct(outputs[0][0][0])

    npflux = flatten_floatstruct(outputs[2][0][0])
    for val in np.hstack([np.array(npspe[1:]).reshape(7,8), np.array(npflux).reshape(7, 3*4)]):
        species = Species(charge                           =val[0] ,
                          mass                             =val[1] ,
                          density                          =val[2] ,
                          density_log_gradient             =val[3] ,
                          temperature                      =val[4] ,
                          temperature_log_gradient         =val[5] ,
                          toroidal_velocity                =val[6] ,
                          toroidal_velocity_gradient       =val[7] ,
                          particle_phi_potential           =val[8] ,
                          particle_a_parallel              =val[9] ,
                          particle_b_field_parallel        =val[10],
                          heat_phi_potential               =val[11],
                          heat_a_parallel                  =val[12],
                          heat_b_field_parallel            =val[13],
                          kinetic_momentum_phi_potential   =val[14],
                          kinetic_momentum_a_parallel      =val[15],
                          kinetic_momentum_b_field_parallel=val[16],
                          field_momentum_phi_potential     =val[17],
                          field_momentum_a_parallel        =val[18],
                          field_momentum_b_field_parallel  =val[19],
                          point=point)
        species.save()

    code = matdict[6][0][0]
    if np.all(code[3] == -999999999):
        output_flag = None
    elif np.all(code[3] == 1):
        output_flag = True
    else:
        raise NotImplementedError

    try:
        name = code[0][0]
    except IndexError:
        name = None
    try:
        version = code[1][0]
    except IndexError:
        version = None 
    try:
        parameters = code[2][0]
    except IndexError:
        parameters = None
    code = Code(point=point, name=name, version=version, parameters=parameters, output_flag=output_flag)
    code.save()

#purge_tables()
matgkdb = io.loadmat('gkdb.mat')
with open('output.tab', newline='') as f:
    reader = csv.reader(f)
    reader.__next__()
    for ii, (matdict, eigenfunc_line) in enumerate(zip(matgkdb['gyrokinetic_linear'][0], reader)):
        matdict_to_SQL(matdict, eigenfunc_line)
