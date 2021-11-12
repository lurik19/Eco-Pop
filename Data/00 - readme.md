---------- run*.txt files structure ----------
Vectors which depend on the resource, i
row #0: beta_i (vector)
row #1: mean fitness of the individuals that eat a resource (<e^Y>_i) (vector)
row #2: mean LOG fitness of the individuals that eat a resource: <Y>_i (vector)
row #3: F_i (vector)

Vectors which depend on the ecotype, b
row #4: mean fitness of every ecotype (e^Y_b) (vector)
row #5: mean LOG fitness of every ecotype: Y_b (vector)
row #6: genome size of every ecotype (\sum_i a_{b, i}) (vector)
row #7: number of individuals in the ecotypes (vector)
row #8: hashing value of every ecotype (vector)

Scalars
row #9: mean fitness of the entire population (<e^Y>) (scalar)
row #10: mean LOG fitness of the entire population (<Y>) (scalar)
row #11: Coefficient of variation of e^Y over the ecotypes (scalar)
row #12: Coefficient of variation of e^Y over the strains (scalar)
row #13: Coefficient of variation CoV(Y) over the ecotypes (scalar)
row #14: Coefficient of variation CoV(Y) over the strains (scalar)
row #15: fraction of generalist strain, R - <k_\mu> (scalar, scalar)
row #16: standard deviation of Y_i (scalar)
row #17: scaled mutation fitness rate (scalar)
row #18: number of sequenced species (scalar)

The first 19 rows describe the first simulation (with the given parameters), the second 19 rows the second simulation and so on.

---------- n_gen*.txt files structure ----------
if (fast == 0)
	Vectors which depend on the resource, i
	row #0: beta_i (vector)
	row #1: mean fitness of the individuals that eat a resource (<e^Y>_i) (vector)
	row #2: mean LOG fitness of the individuals that eat a resource: <Y>_i (vector)
	row #3: F_i (vector)

	Vectors which depend on the ecotype, b
	row #4: mean fitness of every ecotype (e^Y_b) (vector)
	row #5: mean LOG fitness of every ecotype: Y_b (vector)
	row #6: genome size of every ecotype (\sum_i a_{b, i}) (vector)
	row #7: number of individuals in the ecotypes (vector)
	row #8: hashing value of every ecotype (vector)

	Scalars
	row #9: mean fitness of the entire population (<e^Y>) (scalar)
	row #10: mean LOG fitness of the entire population (<Y>) (scalar)
	row #11: Coefficient of variation of e^Y over the ecotypes (scalar)
	row #12: Coefficient of variation of e^Y over the strains (scalar)
	row #13: Coefficient of variation CoV(Y) over the ecotypes (scalar)
	row #14: Coefficient of variation CoV(Y) over the strains (scalar)
	row #15: fraction of generalist strain, R - <k_\mu> (scalar, scalar)
	row #16: number of sequenced species (scalar)
	row #17: Lyapunov function (scalar)
	
if (fast == 1)
	the file is not created

Every dt_gen we produce 18 rows.
The first 18 * n_gen/dt_gen rows describe the first simulation, the second 18 * n_gen/dt_gen rows the second simulation and so on.
