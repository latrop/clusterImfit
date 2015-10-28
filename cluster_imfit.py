#! /usr/bin/env python

import sys
import os
from os import getcwd
from multiprocessing import Pool
import subprocess
from shutil import move

from numpy import argmin

from libs.read_input import ImfitModel, GeneralParams
from libs.pygene.gene import FloatGene, FloatGeneMax
from libs.pygene.organism import Organism, MendelOrganism
from libs.pygene.population import Population


def remove(pth):
    if os.path.exists(pth):
        os.remove(pth)


def run_imfit_parallel(runString):
    chisq = 1e10
    proc = subprocess.Popen(runString, stdout=subprocess.PIPE, shell=True)
    for line in proc.stdout:
        if "Reduced Chi^2 =" in line:
            chisq = float(line.split()[3])
    proc.wait()
    return chisq


class Converger(MendelOrganism):
    """
    Implements the organism which tries to converge a function
    """
    model = ImfitModel(sys.argv[1])
    genome = model.create_genome()
    def prepare_fitness(self):
        genome = {}
        for key in self.genes.keys():
            genome[key] = self[key]
        self.model.genome_to_model(genome)
        fname = self.model.create_input_file(fixAll=True)
        runString = "%simfit -c %s %s " % (gParams.imfitPath, fname, gParams.fitsToFit)
        runString += " --nm --max-threads 1 "
        runString += " --save-params /dev/null "
        if gParams.PSF != "none":
            runString += " --psf %s " % (gParams.PSF)
        if gParams.mask != "none":
            runString += " --mask %s " % (gParams.mask)
        if gParams.weight != "none":
            runString += " --noise %s " % (gParams.weight)
        if gParams.readNoise != "none":
            runString += " --readnoise=%1.2f " % (gParams.readNoise)
        if gParams.gain != "none":
            runString += " --gain=%1.2f " % (gParams.readNoise)
        runString += " %s " % (gParams.addImfitStr)
        runString += " ; rm %s " % (fname)
        result = pool.apply_async(run_imfit_parallel, [runString])
        self.chisq = result

    def fitness(self):
        return self.chisq.get()

    def save_results(self, outFile):
        fname = self.model.create_input_file(fixAll=True)
        runString = "%smakeimage %s --refimage %s " % (gParams.imfitPath, fname, gParams.fitsToFit)
        if gParams.PSF != "none":
            runString += " --psf %s " % (gParams.PSF)
        runString += "--output %s ; rm %s" % (outFile, fname)
        proc = subprocess.Popen(runString, stdout=subprocess.PIPE, shell=True)
        proc.wait()

    def run_lm_optimisation(self, ident):
        fname = "%s/results/%i_lm_input.dat" % (getcwd(), ident)
        genome = {}
        for key in self.genes.keys():
            genome[key] = self[key]
        self.model.genome_to_model(genome)
        self.model.create_input_file(fname)
        runString = "%simfit -c %s %s " % (gParams.imfitPath, fname, gParams.fitsToFit)
        runString += " --max-threads 1 "
        if gParams.PSF != "none":
            runString += " --psf %s " % (gParams.PSF)
        if gParams.mask != "none":
            runString += " --mask %s " % (gParams.mask)
        if gParams.weight != "none":
            runString += " --noise %s " % (gParams.weight)
        if gParams.readNoise != "none":
            runString += " --readnoise=%1.2f " % (gParams.readNoise)
        if gParams.gain != "none":
            runString += " --gain=%1.2f " % (gParams.readNoise)
        runString += "--ftol 0.000001"
        runString += " --save-params %s/results/%i_lm_result.dat " % (getcwd(), ident)
        runString += " --save-model %s/results/%i_lm_model.fits " % (getcwd(), ident)
        runString += " --save-residual %s/results/%i_lm_residual.fits " % (getcwd(), ident) 
        runString += " %s " % (gParams.addImfitStr)
        # runString += " ; rm %s" % (fname)
        result = pool.apply_async(run_imfit_parallel, [runString])
        return result


if __name__ == '__main__':
    gParams = GeneralParams("config.dat")
    pop = Population(species=Converger, init=gParams.zeroGenSize,
                     childCount=gParams.popSize,
                     childCull=gParams.selectNbest,
                     numNewOrganisms=gParams.addNew)

    if not os.path.exists("%s/results/generations/" % getcwd()):
        os.makedirs("%s/results/generations/" % getcwd())
    pool = Pool(gParams.numOfCores)

    iGen = 0
    bestFitness = []
    avgFitness = []
    print "Starting genetic algorithm"
    while 1:
        best = pop.best()
        ftns = best.get_fitness()
        avgFtns = pop.fitness()
        print "generation %i: best=%5.2f average=%5.2f" % (iGen, ftns, avgFtns),
        bestFitness.append(ftns)
        avgFitness.append(avgFtns)
        best.save_results("%s/results/generations/gen_%03i.fits" % (getcwd(), iGen))
        if (iGen > gParams.fSpan):
            relBestFitnessChange = abs(bestFitness[-1] - bestFitness[-gParams.fSpan]) / bestFitness[-1]
            relAvgFitnessChange = abs(avgFitness[-1]-avgFitness[-gParams.fSpan]) / avgFitness[-1]
            print " (delta=%1.3e)" % (max(relBestFitnessChange, relAvgFitnessChange))
            if  (relBestFitnessChange < gParams.fTol) and (relAvgFitnessChange < gParams.fTol):
                print "\n Method converged"
                break
        else:
            print
        if iGen >= gParams.maxGenNumber:
            print "\n Maximum number of generation reached."
            break
        iGen += 1
        pop.gen()

    bestOrganisms = sorted(pop)[0:gParams.numOfCores]
    print "\n Starting L-M optimisation"
    result = [org.run_lm_optimisation(i) for i,org in enumerate(bestOrganisms)]
    chiSqValues = [r.get() for r in result]
    bestModelNumber = argmin(chiSqValues)
    print "Chi.sq. of the best model = %1.3f" % (chiSqValues[bestModelNumber])
    # Remove all models except the best one
    for i in xrange(gParams.numOfCores):
        if i != bestModelNumber:
            remove("%s/results/%i_lm_input.dat" % (getcwd(), i))
            remove("%s/results/%i_lm_result.dat" % (getcwd(), i))
            remove("%s/results/%i_lm_model.fits" % (getcwd(), i))
            remove("%s/results/%i_lm_residual.fits" % (getcwd(), i))
        else:
            move("%s/results/%i_lm_input.dat" % (getcwd(), i), "%s/results/input.dat" % (getcwd()))
            move("%s/results/%i_lm_result.dat" % (getcwd(), i), "%s/results/result.dat" % getcwd())
            move("%s/results/%i_lm_model.fits" % (getcwd(), i), "%s/results/model.fits" % getcwd())
            move("%s/results/%i_lm_residual.fits" % (getcwd(), i), "%s/results/residual.fits" % getcwd())
