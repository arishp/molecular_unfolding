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
hubo = {('x_0_0', 'x_0_1'): 40.0, ('x_0_0', 'x_0_2', 'x_1_0', 'x_1_2'): -1.92, ('x_0_0', 'x_0_2', 'x_1_0'): 7.68, ('x_0_0', 'x_0_2', 'x_1_1', 'x_1_3'): -9.6, 
('x_0_0', 'x_0_2', 'x_1_2'): -7.68, ('x_0_0', 'x_0_2'): 59.36, ('x_0_0', 'x_0_3'): 40.0, ('x_0_0', 'x_1_0'): 0.8, ('x_0_0', 'x_1_2'): -0.8, ('x_0_0',): -34.8, 
('x_0_1', 'x_0_1'): 10.32, ('x_0_1', 'x_0_2'): 40.0, ('x_0_1', 'x_0_3', 'x_1_0', 'x_1_2'): -1.92, ('x_0_1', 'x_0_3', 'x_1_0'): 7.68, ('x_0_1', 'x_0_3', 'x_1_1', 'x_1_3'): -9.6, 
('x_0_1', 'x_0_3', 'x_1_2'): -7.68, ('x_0_1', 'x_0_3'): 59.36, ('x_0_1', 'x_1_1'): -1.78885438199983, ('x_0_1', 'x_1_3'): 1.78885438199983, ('x_0_1',): -40.0, 
('x_0_2', 'x_0_2'): 10.32, ('x_0_2', 'x_0_3'): 40.0, ('x_0_2', 'x_1_0'): -0.8, ('x_0_2', 'x_1_2'): 0.8, ('x_0_2',): -45.2, ('x_0_3', 'x_0_3'): 10.32, 
('x_0_3', 'x_1_1'): 1.78885438199983, ('x_0_3', 'x_1_3'): -1.78885438199983, ('x_0_3',): -40.0, ('x_1_0', 'x_1_0'): 18.08, ('x_1_0', 'x_1_1'): 40.0, ('x_1_0', 'x_1_2'): 43.84, 
('x_1_0', 'x_1_3'): 40.0, ('x_1_0',): -32.16, ('x_1_1', 'x_1_1'): 20.0, ('x_1_1', 'x_1_2'): 40.0, ('x_1_1', 'x_1_3'): 40.0, ('x_1_1',): -40.0, ('x_1_2', 'x_1_2'): 18.08, 
('x_1_2', 'x_1_3'): 40.0, ('x_1_2',): -47.84, ('x_1_3', 'x_1_3'): 20.0, ('x_1_3',): -40.0, (): 22.58}

bqm = dimod.make_quadratic(hubo, 12.0, dimod.BINARY)

sampler = neal.SimulatedAnnealingSampler()
sample_size=10
sampleset = sampler.sample(bqm, num_reads=sample_size)
sa_solution = sampleset.first.sample
print("\nBEST SA RESULT:\n---- -- ------")
items = sa_solution.keys()
for item in items:
    if sa_solution[item] == 1:
        print(item, ':', sa_solution[item])

sampler = EmbeddingComposite(DWaveSampler())
sampleset = sampler.sample(bqm, num_reads=1000)
qa_solution = sampleset.first.sample
print("\nBEST QA RESULT:\n---- -- ------")
items = qa_solution.keys()
for item in items:
    if qa_solution[item] == 1:
        print(item, ':', qa_solution[item])
