# hubo[('x00', 'x00')] = 18.0
# hubo[('x00', 'x01')] = 40.0
# hubo[('x00', 'x02')] = 44.0
# hubo[('x00', 'x03')] = 40.0
# hubo[('x00', )] = -38.0
# hubo[('x01', 'x01')] = 18.0
# hubo[('x01', 'x02')] = 40.0
# hubo[('x01', 'x03')] = 44.0
# hubo[('x01', )] = -40.0
# hubo[('x02', 'x02')] = 18.0
# hubo[('x02', 'x03')] = 40.0
# hubo[('x02', )] = -42.0
# hubo[('x03', 'x03')] = 18.0
# hubo[('x03', )] = -40.0
# hubo[()] =  12.75
# print("\nHUBO:\n", hubo)

import dimod
import neal
from dwave.system import DWaveSampler, EmbeddingComposite
import sympy as sp

thetas = {0: 0, 1: sp.pi/2, 2: sp.pi, 3: 3*sp.pi/2}
sp.pprint(thetas)
A_const = sp.Symbol('A_const')
A_const = 100.0
hubo = {
    ('x_0_0', 'x_0_1'): 2*A_const,
    ('x_0_0', 'x_0_2'): 2*A_const + 19.36,
    ('x_0_0', 'x_0_3'): 2*A_const,
    ('x_0_0', 'x_1_0'): 0.8,
    ('x_0_0', 'x_1_1'): 0,
    ('x_0_0', 'x_1_2'): -0.8,
    ('x_0_0', 'x_1_3'): 0,
    ('x_0_0',): -2*A_const + 5.2, 
    ('x_0_1', 'x_0_2'): 2*A_const, 
    ('x_0_1', 'x_0_3'): 2*A_const + 19.36, 
    ('x_0_1', 'x_1_0'): 0,
    ('x_0_1', 'x_1_1'): -1.78885438199983, 
    ('x_0_1', 'x_1_2'): 0, 
    ('x_0_1', 'x_1_3'): 1.78885438199983, 
    ('x_0_1',): 2*A_const, 
    ('x_0_2', 'x_0_3'): 2*A_const, 
    ('x_0_2', 'x_1_0'): -0.8, 
    ('x_0_2', 'x_1_1'): 0, 
    ('x_0_2', 'x_1_2'): 0.8, 
    ('x_0_2', 'x_1_3'): 0, 
    ('x_0_2',): -2*A_const - 5.2, 
    ('x_0_3', 'x_1_0'): 0, 
    ('x_0_3', 'x_1_1'): 1.78885438199983, 
    ('x_0_3', 'x_1_2'): 0,
    ('x_0_3', 'x_1_3'): -1.78885438199983, 
    ('x_0_3',): - 2*A_const, 
    ('x_1_0', 'x_1_1'): 2*A_const, 
    ('x_1_0', 'x_1_2'): 2*A_const + 3.84, 
    ('x_1_0', 'x_1_3'): 2*A_const, 
    ('x_1_0',): -2*A_const + 7.84, 
    ('x_1_1', 'x_1_2'): 2*A_const, 
    ('x_1_1', 'x_1_3'): 2*A_const, 
    ('x_1_1',): -2*A_const,
    ('x_1_2', 'x_1_3'): 2*A_const, 
    ('x_1_2',): -2*A_const - 7.84,  
    ('x_1_3',): -2*A_const, 
    (): 2*A_const - 17.42}

sp.pprint(hubo)

bqm = dimod.make_quadratic(hubo, 12.0, dimod.BINARY)

sampler = neal.SimulatedAnnealingSampler()
sample_size=1000
sampleset = sampler.sample(bqm, num_reads=sample_size)
sa_solution = sampleset.first.sample
print(sa_solution)
print("\nBEST SA RESULT:\n---- -- ------")
items = sa_solution.keys()
for item in items:
    if sa_solution[item] == 1:
        print(item, ':', sa_solution[item])

sampler = EmbeddingComposite(DWaveSampler())
sampleset = sampler.sample(bqm, num_reads=1000)
qa_solution = sampleset.first.sample
print(qa_solution)
print("\nBEST QA RESULT:\n---- -- ------")
items = qa_solution.keys()
for item in items:
    if qa_solution[item] == 1:
        print(item, ':', qa_solution[item])
