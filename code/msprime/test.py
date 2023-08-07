import tskit
import msprime
import stdpopsim

ts = msprime.sim_ancestry(5, model = "hudson", random_seed= 42)

print(ts.draw_text(y_axis=True, y_ticks=[0, 1, 2, 3, 4, 5]))

ts1 = msprime.sim_ancestry(2, sequence_length=10, recombination_rate=0.1, random_seed=42)
print(ts1.draw_text(y_axis=True, y_ticks=[0, 1, 2, 3, 4, 5]))

