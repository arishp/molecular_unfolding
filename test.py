import dimod
import neal
from dwave.system import DWaveSampler, EmbeddingComposite


hubo = {}

hubo[('x00', 'x00')] = 18.0
hubo[('x00', 'x01')] = 40.0
hubo[('x00', 'x02')] = 44.0
hubo[('x00', 'x03')] = 40.0
hubo[('x00', )] = -38.0
hubo[('x01', 'x01')] = 18.0
hubo[('x01', 'x02')] = 40.0
hubo[('x01', 'x03')] = 44.0
hubo[('x01', )] = -40.0
hubo[('x02', 'x02')] = 18.0
hubo[('x02', 'x03')] = 40.0
hubo[('x02', )] = -42.0
hubo[('x03', 'x03')] = 18.0
hubo[('x03', )] = -40.0
hubo[()] =  12.75
print("\nHUBO:\n", hubo)

bqm = dimod.make_quadratic(hubo, 12.0, dimod.BINARY)

sampler = neal.SimulatedAnnealingSampler()
sample_size=10
sampleset = sampler.sample(bqm, num_reads=sample_size)
sa_solution = sampleset.first.sample
print("\nSA RESULTS:\n",sampleset)
print("\nBEST SA RESULT:\n---- -- ------\n", sa_solution)

sampler = EmbeddingComposite(DWaveSampler())
sampleset = sampler.sample(bqm, num_reads=1000)
qa_solution = sampleset.first.sample
print("\nQA RESULTS:\n",sampleset)
print("\nBEST QA RESULT:\n---- -- ------\n", qa_solution)
