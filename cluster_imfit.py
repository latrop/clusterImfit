#! /usr/bin/env python

import copy

import os
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
    runString = "./imfit -c %s %s " % (fname, gParams.fitsToFit)
    runString += "--readnoise=%1.2f " % (gParams.readNoise)
    runString += " --nm --max-threads 1 "
    runString += " --save-params /dev/null "
    if gParams.PSF != "none":
        runString += " --psf %s " % (gParams.PSF)
    if gParams.mask != "none":
        runString += " --mask %s" % (gParams.mask)
    runString += " ; rm %s " % (fname)
    proc = subprocess.Popen(runString, stdout=subprocess.PIPE, shell=True)
    for line in proc.stdout:
        if "Reduced Chi^2 =" in line:
            chisq = float(line.split()[3])
    proc.wait()
    return chisq


#def create_Converger_class(model):
class Converger(MendelOrganism):
    """
    Implements the organism which tries to converge a function
    """
    model = ImfitModel("input.dat")
    genome = model.create_genome()
    def prepare_fitness(self):
        genome = {}
        for key in self.genes.keys():
            genome[key] = self[key]
        self.model.genome_to_model(genome)
        self.result = pool.apply_async(calculate_fitness, [self.model, gParams])

    def fitness(self):
        return self.result.get()

    def save_results(self, outFile):
        fname = self.model.create_input_file(fixAll=True)
        runString = "./makeimage %s --refimage %s " % (fname, gParams.fitsToFit)
        if gParams.PSF != "none":
            runString += " --psf %s" % (gParams.PSF)
        runString += "--output %s ; rm %s" % (outFile, fname)
        proc = subprocess.Popen(runString, stdout=subprocess.PIPE, shell=True)
        proc.wait()

    def run_lm_optimisation(self, ident):
        fname = "./results/%i_lm_input.dat" % (ident)
        self.model.create_input_file(fname)
        runString = "./imfit -c %s %s " % (fname, gParams.fitsToFit)
        runString += "--readnoise=%1.2f " % (gParams.readNoise)
        # runString += "--max-threads 1 "
        if gParams.PSF != "none":
            runString += " --psf %s " % (gParams.PSF)
        if gParams.mask != "none":
            runString += " --mask %s " % (gParams.mask)
        runString += "--ftol 0.001"
        runString += " --save-params ./results/%i_lm_result.dat " % (ident)
        runString += " --save-model ./results/%i_lm_model.fits " % (ident)
        runString += " --save-residual ./results/%i_lm_residual.fits " % (ident) 
        # runString += " ; rm %s" % (fname)
        proc = subprocess.Popen(runString, stdout=subprocess.PIPE, shell=True)
        for line in proc.stdout:
            print line


def main():
    iGen = 0
    while iGen<=gParams.maxGenNumber:
        best = pop.best()
        print "generation %s: %s best=%s average=%s" % (
            iGen, repr(best), best.get_fitness(), pop.fitness())
        best.save_results("./results/generations/gen_%03i.fits" % (iGen))
        if best.get_fitness() <= 0:
            print "cracked!"
            break
        iGen += 1
        pop.gen()
    best.run_lm_optimisation(iGen)

if __name__ == '__main__':
    gParams = GeneralParams("config.dat")
    pop = Population(species=Converger, init=gParams.zeroGenSize,
                     childCount=gParams.popSize,
                     childCull=gParams.selectNbest,
                     numNewOrganisms=gParams.addNew)

    if not os.path.exists("./results/generations/"):
        os.makedirs("./results/generations/")
    pool = Pool(gParams.numOfCores)
    main()
