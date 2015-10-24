#! /usr/bin/env python

from os import remove
from multiprocessing import Pool
import subprocess

from pygene.gene import FloatGene, FloatGeneMax
from pygene.organism import Organism, MendelOrganism
from pygene.population import Population

from libs.read_input import ImfitModel, GeneralParams

def calculate_fitness(model, gParams):
    """
    Implements the 'fitness function' for this species.
    Organisms try to evolve to minimise this function's value
    """
    fname = model.create_input_file(fixAll=True)
    runString = "./imfit -c %s %s --readnoise=4.3 --nm --max-threads 1 ; rm %s" % (fname, gParams.fitsToFit, fname)
    proc = subprocess.Popen(runString, stdout=subprocess.PIPE, shell=True)
    for line in proc.stdout:
        if "CHI-SQUARE =" in line:
            chisq = float(line.split()[2])
    proc.wait()
    return chisq


#def create_Converger_class(model):
class Converger(MendelOrganism):
    """
    Implements the organism which tries to converge a function
    """
    model = ImfitModel("input.dat")
    gParams = GeneralParams("config.dat")
    genome = model.create_genome()
    def prepare_fitness(self):
        genome = {}
        for key in self.genes.keys():
            genome[key] = self[key]
        self.model.genome_to_model(genome)
        self.result = pool.apply_async(calculate_fitness, [self.model, self.gParams])

    def fitness(self):
        return self.result.get()


class GalaxyPopulation(Population):

    initPopulation = 10
    species = Converger
    
    # cull to this many children after each generation
    childCull = 10
    
    # number of children to create after each generation
    childCount = 50
    
    mutants = 0.25


#model = ImfitModel("input.dat")
#gParams = GeneralParams("config.dat")
ph = GalaxyPopulation()
def main(nfittest=10, nkids=100):
    i = 0
    while i<=0:
        b = ph.best()
        print "generation %s: %s best=%s average=%s)" % (
            i, repr(b), b.get_fitness(), ph.fitness())
        if b.get_fitness() <= 0:
            print "cracked!"
            break
        i += 1
        ph.gen()



if __name__ == '__main__':
    pool = Pool(3)
    main()
