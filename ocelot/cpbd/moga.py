import random
from deap import base
from deap import creator
from deap import tools

import pickle
import numpy as np
from scipy.optimize import *

try:
    from mpi4py import MPI

    MPI_COMM = MPI.COMM_WORLD
    MPI_SIZE = MPI_COMM.Get_size()
    MPI_RANK = MPI_COMM.Get_rank()

except Exception:
    MPI_SIZE = 1
    MPI_RANK = 0


class Moga():

    def __init__(self, bounds, weights=(-1.0, -1.0)):

        # population and elite population sizes
        self.n_pop = 100
        self.elite_num = int(self.n_pop / 10)
        self.new_num = int(self.n_pop / 10)

        # number of generation
        self.n_gen = 10

        # infinite value
        self.inf_val = float("inf")
        self.seed = None

        self.c_iter = None
        self.c_ind = None

        self.weights = weights
        self.problem_size = len(self.weights)

        self.penalty = None

        self.log_print = True if MPI_RANK == 0 else False
        self.log_file = 'moga_result.dat' if MPI_RANK == 0 else None
        self.plt_file = 'moga_plot.dat' if MPI_RANK == 0 else None

        self.fit_func = lambda x: None
        self.fit_func_args = []

        self.vars_num = len(bounds)
        self.bounds_min = []
        self.bounds_max = []
        for i in range(len(bounds)):
            self.bounds_min.append(bounds[i][0])
            self.bounds_max.append(bounds[i][1])

    def set_params(self, n_pop=None, weights=None, elite=None, penalty=None, n_gen=None, seed=None, log_print=None,
                   log_file=False, plt_file=False):

        if n_pop != None:
            self.n_pop = n_pop
            self.elite_num = int(self.n_pop / 10)
            self.new_num = int(self.n_pop / 10)

        if weights != None:
            self.weights = weights
            self.problem_size = len(self.weights)

        if elite != None and elite < self.n_pop:
            self.elite = elite

        if penalty != None:
            self.penalty = penalty

        if n_gen != None:
            self.n_gen = n_gen

        if seed != None:
            self.seed = seed

        if log_print != None and MPI_RANK == 0:
            self.log_print = True

        if log_file != False and MPI_RANK == 0:
            self.log_file = log_file

        if plt_file != False and MPI_RANK == 0:
            self.plt_file = plt_file

    def generate_ind(self):
        return [np.random.uniform(self.bounds_min[i], self.bounds_max[i]) for i in range(self.vars_num)]

    def init_deap_functions(self):

        creator.create("Fitness", base.Fitness, weights=self.weights)
        creator.create("Individual", list, fitness=creator.Fitness)

        self.toolbox = base.Toolbox()

        self.toolbox.register("individual", tools.initIterate, creator.Individual, self.generate_ind)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)

        self.toolbox.register("evaluate", self.fit_func)

        if self.penalty != None:
            self.toolbox.decorate("evaluate", tools.DeltaPenality(self.feasible, self.inf_val))

    def feasible(self):
        return True

    def optimize(self, init_pop):

        # Create initial population
        pop = self.init_pop(init_pop)

        if self.log_print: print("-- Iteration 0 --")
        if self.log_file:
            fh1 = open(self.log_file, 'w')
            fh1.close()

        # Evaluate initial population
        self.c_iter = 0
        pop = self.eval_pop(pop)

        # This is just to assign the crowding distance to the individuals (no actual selection is done)
        pop = self.toolbox.select(pop, len(pop))

        # Begin the evolution
        for g in range(self.n_gen):

            # Select good and bad individuals
            good_inds = self.get_good_inds(pop)

            # Create elite population (non dominated individuals)
            nond_inds = self.get_nondominated_inds(pop)

            self.c_iter = g + 1
            if self.log_print: print("-- Iteration %i --" % self.c_iter)

            # Save current population to file
            if self.plt_file != None:

                data_file = []
                data_file.append([g + 1, self.n_gen])

                val_g = []
                for ind in good_inds:
                    val_g.append(ind.fitness.values)
                data_file.append(val_g)

                val_nd = []
                for ind in nond_inds:
                    val_nd.append(ind.fitness.values)
                data_file.append(val_nd)

                with open(self.plt_file, 'wb') as fh2:
                    pickle.dump(data_file, fh2, protocol=2)

            if self.log_file != None and (len(good_inds) > 0):
                fh1 = open(self.log_file, 'a')
                fh1.write("\nBest individuals\n")
                for ind in pop:
                    fh1.write("ind --> fit_func: " + str(ind) + ' --> ' + str(ind.fitness.values) + '\n')
                fh1.close()

            # Vary the population
            offspring = self.apply_preselect(pop)

            # Apply crossover on the offspring
            offspring = self.apply_crossover(offspring)

            # Apply mutation on the offspring
            offspring = self.apply_mutation(offspring)

            # Select individuals with invalid fitness and evaluate
            if MPI_RANK == 0:
                invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            else:
                invalid_ind = None

            invalid_ind = self.eval_pop(invalid_ind)

            # Select the next generation population
            pop = self.apply_select(pop, offspring)

            # Replace random or bad individuals by new random and evaluate
            if self.n_pop - len(good_inds) < self.new_num:
                pop = self.replace_rand_by_new_inds(pop, self.new_num)
            else:
                pop = self.replace_bad_by_rand_inds(pop)

            if MPI_RANK == 0:
                invalid_ind = [ind for ind in pop if not ind.fitness.valid]
            else:
                invalid_ind = None

            invalid_ind = self.eval_pop(invalid_ind)

        self.c_iter = None

        return pop

    def init_pop(self, init_pop):

        if MPI_RANK != 0: return None

        pop = self.toolbox.population(n=self.n_pop)
        if init_pop != None:

            init_pop_len = len(init_pop) if len(init_pop) <= self.n_pop else self.n_pop

            for i in range(init_pop_len):
                for j in range(len(init_pop[i])):
                    pop[i][j] = init_pop[i][j]

        return pop

    def eval_pop(self, pop):

        # Split data
        if MPI_SIZE > 1:
            data_in = None
            keys = None
            if MPI_RANK == 0:
                data_in = []
                keys = []
                length = int(len(pop) / MPI_SIZE)
                for i in range(MPI_SIZE):
                    start = i * length
                    stop = start + length
                    data_in.append(pop[start:stop])
                    keys.append([ii for ii in range(start, stop)])

                j = 0
                for i in range(stop, len(pop)):
                    data_in[j] += [pop[i]]
                    keys[j] += [i]
                    j += 1

            data_in = MPI_COMM.scatter(data_in, root=0)
            keys = MPI_COMM.scatter(keys, root=0)
        else:
            data_in = pop
            keys = [ii for ii in range(len(data_in))]

        # Evaluate initial population
        fitnesses = [self.toolbox.evaluate(x0=x, iter_data=(key, self.c_iter), args=self.fit_func_args) for key, x in
                     zip(keys, data_in)]

        # Merge data
        if MPI_SIZE > 1:
            data_out = MPI_COMM.gather(fitnesses, root=0)

            if MPI_RANK == 0:
                fitnesses = []
                for i in data_out:
                    fitnesses += i[:length]
                for i in data_out:
                    fitnesses += i[length:length + 1]
                    j -= 1
                    if j == 0: break

        # Update population data
        if MPI_RANK == 0:
            for ind, fit in zip(pop, fitnesses):
                ind.fitness.values = fit
        else:
            pop = None

        if self.log_print: print("Evaluated %i" % (len(pop)))

        return pop

    def get_good_inds(self, pop):

        if MPI_RANK != 0: return None

        g_inds = []
        gb = [0, 0]

        for ppp in pop:

            good_val_checker = True
            for j in range(self.problem_size):
                if ppp.fitness.values[j] == self.inf_val:
                    good_val_checker = False

            if good_val_checker:
                g_inds.append(self.toolbox.clone(ppp))
                gb[0] += 1
            else:
                gb[1] += 1

        if self.log_print:
            print("good/bad solitions %i / %i" % (gb[0], gb[1]))

        return g_inds

    def get_nondominated_inds(self, pop):

        if MPI_RANK != 0: return None

        nd_inds = tools.sortNondominated(pop, k=len(pop), first_front_only=True)[0]
        nd_inds = [self.toolbox.clone(x) for x in nd_inds]

        if self.log_print:
            print("non dominated %i" % len(nd_inds))

        return nd_inds

    def get_best_inds(self, pop):
        # only for minimization problem

        if MPI_RANK != 0: return None

        best_ind = []
        for j in range(self.problem_size):
            best_ind.append(pop[0])

        for ind in pop:
            for j in range(self.problem_size):
                if ind.fitness.values[j] < best_ind[j].fitness.values[j]:
                    best_ind[j] = ind

        best_ind = [self.toolbox.clone(x) for x in best_ind]

        return best_ind

    def apply_preselect(self, pop):

        if MPI_RANK != 0: return None

        offspring = self.toolbox.preselect(pop, len(pop))
        offspring = [self.toolbox.clone(x) for x in offspring]

        return offspring

    def apply_select(self, pop, offspring):

        if MPI_RANK != 0: return None

        unique = []
        for ooo in offspring:
            if ooo not in pop:
                unique.append(self.toolbox.clone(ooo))

        result = self.toolbox.select(pop + unique, self.n_pop)

        return result

    def apply_crossover(self, pop):

        if MPI_RANK != 0: return None

        for child1, child2 in zip(pop[::2], pop[1::2]):

            # cross two individuals with probability CXPB
            if np.random.random() < self.cxpb:
                self.toolbox.mate(child1, child2)

                # fitness values of the children must be recalculated later
                del child1.fitness.values
                del child2.fitness.values

        return pop

    def apply_mutation(self, pop):

        if MPI_RANK != 0: return None

        for mutant in pop:

            # mutate an individual with probability MUTPB
            if np.random.random() < self.mutpb:
                self.toolbox.mutate(mutant)
                del mutant.fitness.values

        return pop

    def replace_rand_by_good_inds(self, pop, good_inds, num=None):

        if MPI_RANK != 0: return None

        n_good = len(good_inds)

        if num == None:
            num = n_good
        else:
            num = np.min((num, n_good))

        i_pop = range(num)
        i_good = range(num)

        for i_p, i_g in zip(i_pop, i_good):

            if good_inds[i_g] not in pop:
                pop[i_p] = self.toolbox.clone(good_inds[i_g])

        return pop

    def replace_bad_by_rand_inds(self, pop):

        if MPI_RANK != 0: return None

        for ppp in pop:

            good_val_checker = True
            for j in range(self.problem_size):
                if ppp.fitness.values[j] == self.inf_val:
                    good_val_checker = False

            if good_val_checker:
                continue

            new_ind = self.toolbox.individual()

            for j in range(len(ppp)):
                ppp[j] = new_ind[j]

            del ppp.fitness.values

        return pop

    def replace_bad_by_good_inds(self, pop, good_pop):

        if MPI_RANK != 0: return None

        for iii in range(len(pop)):

            if len(good_pop) < 1:
                break

            good_val_checker = True
            for j in range(self.problem_size):
                if pop[iii].fitness.values[j] == self.inf_val:
                    good_val_checker = False

            if good_val_checker:
                continue

            i_nd = len(good_pop) - 1
            if good_pop[i_nd] not in pop:
                pop[iii] = self.toolbox.clone(good_pop[i_nd])
                del good_pop[i_nd]
            else:
                del good_pop[i_nd]

        return pop

    def replace_rand_by_new_inds(self, pop, num):

        if MPI_RANK != 0: return None

        r_pop = tools.selRandom(pop, num)

        nd = self.get_nondominated_inds(pop)
        nd_num = len(nd)
        n_pop_f = int(self.n_pop * 0.3)

        for ppp in r_pop:

            new_ind = self.toolbox.individual()

            if nd_num < n_pop_f and ppp in nd:
                continue

            for j in range(len(ppp)):
                ppp[j] = new_ind[j]

            del ppp.fitness.values

        return pop

    def nsga2(self, fit_func, fit_func_args=[], init_pop=None, cxpb=0.95, mutpb=0.95):

        # init
        random.seed(self.seed)
        np.random.seed(self.seed)

        self.fit_func = fit_func
        self.fit_func_args = fit_func_args
        self.cxpb = cxpb  # cross probability
        self.mutpb = mutpb  # mutation probability

        self.init_deap_functions()
        self.toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=self.bounds_min, up=self.bounds_max,
                              eta=1.0)  # crossover = mate
        self.toolbox.register("mutate", tools.mutPolynomialBounded, eta=10.0, low=self.bounds_min, up=self.bounds_max,
                              indpb=self.mutpb)
        # self.toolbox.register("preselect", tools.selTournamentDCD)
        self.toolbox.register("preselect", tools.selRandom)
        self.toolbox.register("select", tools.selNSGA2)

        if self.log_print: print("Number of used CPU: %i" % MPI_SIZE)

        # optimization
        result = self.optimize(init_pop)

        if self.log_print: print("End of (successful) evolution")

        result_nd = self.get_nondominated_inds(result)

        if self.log_file:
            fh1 = open(self.log_file, 'a')
            fh1.write("\n-------------------------- End of (successful) evolution --------------------------\n")
            fh1.write("\nNon dominated individuals\n")
            for ind in result_nd:
                fh1.write("ind --> fit_func: " + str(ind) + ' --> ' + str(ind.fitness.values) + '\n')
            fh1.close()

        return result_nd
