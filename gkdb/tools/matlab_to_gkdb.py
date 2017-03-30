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

    inputs = matdict[4][0][0]

    collisions = inputs[3][0][0]
    collisions = flatten_floatstruct(collisions)
    collisionality = collisions[0]

    em_effects = inputs[2][0][0]
    em_effects = flatten_floatstruct(em_effects)
    beta = em_effects[0]

    if creator is None:
        return
    point = Point(comment=matdict[0], creator=creator, date=date, beta=beta, collisionality=collisionality)
    point.save()



    equi = inputs[0][0][0]
    sha = equi[4][0][0]
    equi = flatten_floatstruct([equi[0], equi[1], equi[2], equi[3], equi[5], equi[6]])
    npsha = flatten_floatstruct(sha)

    flux_surface = Flux_Surface(point=point,
                                r_minor=equi[0],
                                q=equi[1],
                                magnetic_shear=equi[2],
                                beta_gradient=equi[3],
                                ip_sign=equi[4],
                                b_field_tor_sign=equi[5],
                                c=[npsha[0]] + npsha[2::4],
                                s=[0] + npsha[3::4],
                                dc_dr_minor=[npsha[1]] + npsha[4::4],
                                ds_dr_minor=[0] + npsha[5::4],
                                elongation=1,
                                triangularity=0,
                                squareness=0)
    flux_surface.save(force_insert=True)

    spe = inputs[1][0][0]
    npspe = flatten_floatstruct(spe)

    inp_mode = inputs[4][0][0]
    inp_mode = flatten_floatstruct(inp_mode)
    wavevector = Wavevector(point=point,
                            radial_wavevector=inp_mode[0],
                            binormal_wavevector=inp_mode[1],
                            poloidal_turns=0)
    wavevector.save()
    outputs = matdict[5][0][0]
    eigenvalues = flatten_floatstruct(outputs[0][0][0])
    eigenfunc = flatten_floatstruct(eigenfunc_line[2:])
    if eigenfunc[1] == -9e40:
        poloidal_angle = None
    else:
        poloidal_angle = float(eigenfunc[1])


    for ii, val in enumerate(np.hstack([np.array(eigenvalues).reshape(3,2), np.array(eigenfunc).reshape(3,6)])):
        if np.all([x is not None for x in val]):
            if ii > 1:
                raise Exception('More than one mode!')
            eigenvalue = Eigenvalue(wavevector=wavevector,
                                    growth_rate=val[0],
                                    frequency=val[1],
                                    growth_rate_tolerance=1e-7,
                                    phi_amplitude=1,
                                    phi_parity=1,
                                    a_amplitude=0,
                                    a_parity=1,
                                    b_amplitude=0,
                                    b_parity=1)
                                  
            eigenvalue.save()
            eigenvector = Eigenvector(eigenvalue=eigenvalue,
                                      poloidal_angle=[poloidal_angle],
                                      r_phi_potential_perturbed   =[val[2]], 
                                      i_phi_potential_perturbed   =[val[3]],
                                      r_a_parallel_perturbed      =[val[4]],
                                      i_a_parallel_perturbed      =[val[5]],
                                      r_b_field_parallel_perturbed=[val[6]],
                                      i_b_field_parallel_perturbed=[val[7]])
            eigenvector.save(force_insert=True)
    eigenvalues = flatten_floatstruct(outputs[0][0][0])

    npflux = flatten_floatstruct(outputs[2][0][0])
    for ii, val in enumerate(np.hstack([np.array(npspe[1:]).reshape(7,8), np.array(npflux).reshape(7, 3*4)])):
        if not np.all([x is None for x in val]):
            if ii > 1:
                raise Exception('More than one mode!')
            species = Species(charge                           =val[0] ,
                              mass                             =val[1] ,
                              density                          =val[2] ,
                              temperature                      =val[4] ,
                              toroidal_velocity                =val[6] ,
                              density_log_gradient             =val[3] ,
                              temperature_log_gradient         =val[5] ,
                              toroidal_velocity_gradient       =val[7] ,
                              point=point)
            species.save()

            fluxes = Particle_Fluxes(phi_potential           =val[8] ,
                              a_parallel              =val[9] ,
                              b_field_parallel        =val[10],
                              species=species,
                              eigenvalue=eigenvalue)
            fluxes.save(force_insert=True)
            fluxes = Heat_Fluxes_Lab(
                              phi_potential               =val[11],
                              a_parallel                  =val[12],
                              b_field_parallel            =val[13],
                              species=species,
                              eigenvalue=eigenvalue)
            fluxes.save(force_insert=True)
            fluxes = Momentum_Fluxes_Lab(
                              phi_potential   =val[14],
                              a_parallel      =val[15],
                              b_field_parallel=val[16],
                              species=species,
                              eigenvalue=eigenvalue)
            fluxes.save(force_insert=True)
                              #field_momentum_phi_potential     =val[17],
                              #field_momentum_a_parallel        =val[18],
                              #field_momentum_b_field_parallel  =val[19],

    code = matdict[6][0][0]
    #if np.all(code[3] == -999999999):
    #    output_flag = None
    #elif np.all(code[3] == 1):
    #    output_flag = True
    #else:
    #    raise NotImplementedError

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
    code = Code(point=point,
                name=name,
                version=version,
                parameters=parameters,
                collision_enhancement_factor=collisions[1],
                collision_pitch_only=collisions[2],
                collision_ei_only=collisions[3],
                collision_momentum_conservation=collisions[4],
                collision_energy_conservation=collisions[5],
                collision_finite_larmor_radius=0,
                collision_a_parallel=em_effects[1],
                collision_b_field_parallel=em_effects[2])
    code.save()

#purge_tables()
matgkdb = io.loadmat('gkdb.mat')
with open('output.tab', newline='') as f:
    reader = csv.reader(f)
    reader.__next__()
    for ii, (matdict, eigenfunc_line) in enumerate(zip(matgkdb['gyrokinetic_linear'][0], reader)):
        matdict_to_SQL(matdict, eigenfunc_line)
