#from gkdb.core.model import *
from peewee import PrimaryKeyField, ForeignKeyField, CompositeKey

def get_closest(column, value):
    model_class = column.model_class
    id = model_class.id.alias('id')
    if model_class._meta.primary_key.__class__ in [PrimaryKeyField,
                                                   ForeignKeyField]:
        prim_key_names = [model_class._meta.primary_key.db_column]
    #elif model_class._meta.primary_key.__class__ == CompositeKey:
    #    prim_key_names = [x + '_id' for x in model_class._meta.primary_key.field_names]
    else:
        raise NotImplementedError('Query on table with primary key '
                                  + str(model_class._meta.primary_key.__class__))
    prim_key_fields = [getattr(model_class, name) for name in prim_key_names]
    low_query = (model_class
                 .select(*prim_key_fields, column)
                 .where(column >= value)
                 .order_by(column)
                 .limit(1)
                 )
    high_query = (model_class
                 .select(*prim_key_fields, column)
                 .where(column < value)
                 .order_by(column.desc())
                 .limit(1)
                 )
    join_clause = low_query.union_all(high_query).alias('res')
    
    join_fields = [getattr(join_clause.c, name) for name in prim_key_names]
    query = (model_class
             .select(
                 *join_fields,
                     getattr(join_clause.c, column.name))
             .from_(join_clause)
             .order_by(fn.abs(value - SQL(column.name)))
             .limit(1)
             )
    query_res = query.get()
    res = model_class.select()
    for key in prim_key_fields:
        res = res.where(key == getattr(query_res, key.name))
    return res.get()

#res = get_closest(Wavevector.binormal_wavevector, 0.3)
#print(res.id, res.binormal_wavevector)
#res = get_closest(Flux_Surface.q, 0.3)
#print(res.point_id, res.q)
#res = get_closest(Particle_Fluxes.a_parallel, 0.3)
#print(res.species_id, res.eigenvalue_id, res.a_parallel)
