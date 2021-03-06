#! /usr/bin/env python

import argparse
import datetime
import time
import sys
import os
from os import path
from os import getcwd
from multiprocessing import Pool
import subprocess
from shutil import move
import uuid

from numpy import argmin, isnan

from libs.read_input import ImfitModel, GeneralParams
from libs.pygene.organism import MendelOrganism
from libs.pygene.population import Population


def remove(pth):
    if os.path.exists(pth):
        os.remove(pth)


def run_imfit_parallel(runString, fname=None):
    chisq = 1e10
    stdoutFileName = "%s/results/stdout_%s.dat" % (getcwd(), uuid.uuid4())
    stdoutFile = open(stdoutFileName, "w")
    proc = subprocess.Popen(runString, stdout=stdoutFile, shell=True)
    proc.wait()
    stdoutFile.close()
    for line in open(stdoutFileName):
        if "Reduced Chi^2 =" in line:
            chisq = float(line.split()[3])
    remove(stdoutFileName)
    if fname is not None:
        remove(fname)
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
        imfit_binary = path.join(gParams.imfitPath, "imfit")
        runString = "%s -c %s %s " % (imfit_binary, fname, gParams.fitsToFit)
        runString += " --fitstat-only --max-threads 1 "
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
            runString += " --gain=%1.2f " % (gParams.gain)
        runString += " %s " % (gParams.addImfitStr)
        result = pool.apply_async(run_imfit_parallel, [runString, fname])
        self.chisq = result

    def fitness(self):
        return self.chisq.get()

    def save_results(self, outFile):
        fname = self.model.create_input_file(fixAll=True)
        makeimage_binary = path.join(gParams.imfitPath, "makeimage")
        runString = "%s %s --refimage %s " % (makeimage_binary, fname, gParams.fitsToFit)
        if gParams.PSF != "none":
            runString += " --psf %s " % (gParams.PSF)
        runString += "--output %s" % (outFile)
        proc = subprocess.Popen(runString, stdout=subprocess.PIPE, shell=True)
        proc.wait()
        remove(fname)

    def run_lm_optimisation(self, ident):
        fname = "%s/results/%i_lm_input.dat" % (getcwd(), ident)
        genome = {}
        for key in self.genes.keys():
            genome[key] = self[key]
        self.model.genome_to_model(genome)
        self.model.create_input_file(fname)
        # If we are NOT going to run LM optimisation right now, then
        # save run strings in a file for the future
        imfit_binary = path.join(gParams.imfitPath, "imfit")
        if gParams.runLM == "no":
            script = open("%s/results/run_lm.sh" % (getcwd()), "a")
            runString = "%s -c %i_lm_input.dat ../%s " % (imfit_binary, ident, gParams.fitsToFit)
            if gParams.PSF != "none":
                runString += " --psf ../%s " % (gParams.PSF)
            if gParams.mask != "none":
                runString += " --mask ../%s " % (gParams.mask)
            if gParams.weight != "none":
                runString += " --noise ../%s " % (gParams.weight)
            if gParams.readNoise != "none":
                runString += " --readnoise=%1.2f " % (gParams.readNoise)
            if gParams.gain != "none":
                runString += " --gain=%1.2f " % (gParams.gain)
            runString += "--ftol 0.00001"
            runString += " --save-params %i_lm_result.dat " % (ident)
            runString += " --save-model %i_lm_model.fits " % (ident)
            runString += " --save-residual %i_lm_residual.fits " % (ident)
            runString += " %s \n\n" % (gParams.addImfitStr)
            script.write(runString)
            script.close()
            return
        else:
            runString = "%s -c %s %s " % (imfit_binary, fname, gParams.fitsToFit)
            runString += " --max-threads %i " % (gParams.LMCores)
            if gParams.PSF != "none":
                runString += " --psf %s " % (gParams.PSF)
            if gParams.mask != "none":
                runString += " --mask %s " % (gParams.mask)
            if gParams.weight != "none":
                runString += " --noise %s " % (gParams.weight)
            if gParams.readNoise != "none":
                runString += " --readnoise=%1.2f " % (gParams.readNoise)
            if gParams.gain != "none":
                runString += " --gain=%1.2f " % (gParams.gain)
            runString += "--ftol 0.00001"
            runString += " --save-params %s/results/%i_lm_result.dat " % (getcwd(), ident)
            runString += " --save-model %s/results/%i_lm_model.fits " % (getcwd(), ident)
            runString += " --save-residual %s/results/%i_lm_residual.fits " % (getcwd(), ident)
            runString += " %s \n\n" % (gParams.addImfitStr)
            result = pool.apply_async(run_imfit_parallel, [runString])
            return result


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: ./cluster_imfit.py  model.imfit config.dat")

    logFileName = "%s/log.dat" % (getcwd())
    logFile = open(logFileName, "a", buffering=1)
    logFile.write("\n\n\n################################\n")

    parser = argparse.ArgumentParser()
    parser.add_argument("model")
    parser.add_argument("config")
    for key, value in GeneralParams(sys.argv[2]).params.items():
        parser.add_argument("--%s" % key, default=value, type=type(value))
    gParams = parser.parse_args()
    pop = Population(species=Converger, init=gParams.zeroGenSize,
                     childCount=gParams.popSize,
                     childCull=gParams.selectNbest,
                     numNewOrganisms=gParams.addNew)

    if not os.path.exists("%s/results" % getcwd()):
        os.makedirs("%s/results" % getcwd())
    if (gParams.saveGens == "yes") and (not os.path.exists("%s/results/generations/" % getcwd())):
        os.makedirs("%s/results/generations/" % getcwd())
    pool = Pool(gParams.numOfCores)

    startTime = time.time()
    iGen = 0
    bestFitness = []
    avgFitness = []
    print("Starting genetic algorithm")
    logFile.write("GA optimisation started at %s\n" % datetime.datetime.now().strftime("%d.%m.%Y %H:%M"))
    while 1:
        best = pop.best()
        ftns = best.get_fitness()
        avgFtns = pop.fitness()
        print("generation %i: best=%8.5f average=%8.5f" % (iGen, ftns, avgFtns), end='')
        logFile.write("generation %i: best=%8.5f average=%8.5f" % (iGen, ftns, avgFtns))
        bestFitness.append(ftns)
        avgFitness.append(avgFtns)

        if isnan(ftns) or isnan(avgFtns):
            print("\nNaN values found in model images. Aborting...")
            logFile.write("\nNaN values found in model images. Aborting...")
            exit(1)
        if gParams.saveGens == "yes":
            best.save_results("%s/results/generations/gen_%03i.fits" % (getcwd(), iGen))
        if gParams.genTextFile is not None:
            best.model.model_to_text(iGen, ftns, gParams.genTextFile)
        if (iGen > gParams.fSpan):
            relBestFitnessChange = abs(bestFitness[-1] - bestFitness[-gParams.fSpan]) / bestFitness[-1]
            relAvgFitnessChange = abs(avgFitness[-1]-avgFitness[-gParams.fSpan]) / avgFitness[-1]
            print(" (delta=%1.5e)" % (max(relBestFitnessChange, relAvgFitnessChange)))
            logFile.write(" (delta=%1.5e)\n" % (max(relBestFitnessChange, relAvgFitnessChange)))
            if (relBestFitnessChange < gParams.fTol) and (relAvgFitnessChange < gParams.fTol):
                print("\n GA method converged")
                logFile.write("\n GA method converged\n")
                break
        else:
            print("")
            logFile.write("\n")
        if iGen >= gParams.maxGenNumber:
            print("\n Maximum number of generation reached.")
            logFile.write("\n Maximum number of generation reached.\n")
            break
        iGen += 1
        pop.gen()

    timeSpentSec = time.time() - startTime
    spentTimeString = time.strftime("%Hh:%Mm:%Ss", time.gmtime(timeSpentSec))
    print("Time spent: %s" % spentTimeString)
    logFile.write("Time spent: %s\n" % spentTimeString)

    if gParams.runLM == "no":
        remove("%s/results/run_lm.sh" % getcwd())
    else:
        print("\n Starting L-M optimisation")
        logFile.write("\n Starting L-M optimisation\n")
    bestOrganisms = sorted(pop)[0:gParams.numOfLM]
    result = [org.run_lm_optimisation(i) for i, org in enumerate(bestOrganisms)]
    if gParams.runLM == "yes":
        chiSqValues = [r.get() for r in result]
        bestModelNumber = argmin(chiSqValues)
        print("\nChi.sq. of the best model = %1.3f" % (chiSqValues[bestModelNumber]))
        logFile.write("\nChi.sq. of the best model = %1.3f\n" % (chiSqValues[bestModelNumber]))
        # If LM optimisation was done then remove all models except the best one
        for i in range(gParams.numOfLM):
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

        print("\nChecking boundaries\n")
        # Load resulting model
        resModel = ImfitModel("%s/results/result.dat" % getcwd())
        # Store final parameters to generations file
        if gParams.genTextFile is not None:
            resModel.model_to_text(iGen+1, chiSqValues[bestModelNumber], gParams.genTextFile)
        # Check boundaries
        badParams = best.model.check_boundaries(resModel)
        if badParams:
            fbad = open("%s/results/bad_params.dat" % getcwd(), "w")
            fbad.truncate(0)
            print("Warging: these parameters have values that are close to their boundaries:")
            logFile.write("Warging: some parameters have values that are close to their boundaries")
            logFile.write(" (see 'bad_params.dat').\n")
            for p in badParams:
                print(p)
                fbad.write("%s\n" % p)
        else:
            print("All parameters are inside of their boundaries")
    time.sleep(1)
