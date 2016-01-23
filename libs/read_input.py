#! /usr/bin/env python

import uuid
from os import getcwd
from os.path import exists
import sys
from libs.pygene.gene import FloatGene, FloatGeneMax, IntGene
from libs.pygene.organism import Organism, MendelOrganism


def parse_imfit_line(line):
    """ Function parses line of imfit data file and
    returns parameters"""
    params = line.split()
    name = params[0]
    value = float(params[1]) # Value of the parameter is the second entry of the line
    # (after the name)
    # Lets now find range of values. Its have to contain the coma and must
    # be the second enrty (otherwise imfin wont work).
    # Some parameters can be fixed, so we have to check this possibility at first
    if (len(params) == 2) or ("fixed" in params[2]) :
        # No bounds specified at all or fixed value
        lowerLim = upperLim = None
    else:
        rangeParams = params[2].split(",")
        lowerLim = float(rangeParams[0])
        upperLim = float(rangeParams[1])
	    
    return ImfitParameter(name, value, lowerLim, upperLim)


class ImfitParameter(object):
    """ Just a container of parameter instance:
    parameter name, value and its range"""
    def __init__(self, name, value, lowerLim, upperLim):
        self.name = name
        self.value = value
        self.lowerLim = lowerLim
        self.upperLim = upperLim
        if lowerLim is None:
            self.fixed = True
        else:
            self.fixed = False
    def tostring(self, fixAll):
        if fixAll or self.fixed:
            return "{:12}{:6.2f}     fixed\n".format(self.name, self.value)
        else:
            return "{:12}{:6.2f}{:12.2f},{:1.2f}\n".format(self.name, self.value, self.lowerLim, self.upperLim)
    def change_value(self, newValue):
        self.value = newValue
        if not self.fixed:
            if self.lowerLim > newValue:
                self.lowerLim = newValue - newValue*0.001
            if self.upperLim < newValue:
                self.upperLim = newValue + newValue*0.001


class ImfitFunction(object):
    """ Class represents imfit function with
    all its parameters their ranges"""
    def __init__(self, funcName, ident):
        # ident is a unical number of the function
        # funcName is just a type of the function and it can be 
        # the same for different galaxy components
        self.name = funcName
        self.ident = ident
        self.uname = "%s.%i" % (funcName, ident) # func unique name
        self.params = []
    def add_parameter(self, newParameter):
        self.params.append(newParameter)
    def num_of_params(self):
        return len(self.params)
    def get_par_by_name(self, name):
        for par in self.params:
            if par.name == name:
                return par


class ImfitModel(object):
    """Imfit functions and their parameters"""
    def __init__(self, modelFileName):
        print "Reading '%s':" % (modelFileName)
        # Read imfit input file
        self.listOfFunctions = []
        self.numberOfParams = 0
        funcName = None
        ident = -1
        for line in open(modelFileName):
            sLine = line.strip()
            if sLine.startswith("#"):
                # It is a comment line, just skip it
                continue
            if len(sLine) == 0:
                "Empty line"
                continue
            if "#" in sLine:
                # Drop the comment part of the line if exists
                sLine = sLine[:sLine.index("#")].strip()
            if sLine.startswith("X0"):
                x0 = parse_imfit_line(sLine)
            elif sLine.startswith("Y0"):
                y0 = parse_imfit_line(sLine)
            elif sLine.startswith("FUNCTION"):
                # New function is found.
                ident += 1
                # If we are working already with some function, then
                # the list of parameters for this function is over and we can
                # add it to the function list
                if funcName is not None:
                    self.listOfFunctions.append(currentFunction)
                funcName = sLine.split()[1]
                currentFunction = ImfitFunction(funcName, ident)
                currentFunction.add_parameter(x0)
                currentFunction.add_parameter(y0)
                self.numberOfParams += 2
            else:
                # If line does not contain nor coordinates nor function name
                # then in has to be a parameter line
                param = parse_imfit_line(sLine)
                currentFunction.add_parameter(param)
                self.numberOfParams += 1
        # append the last function
        self.listOfFunctions.append(currentFunction)
        # Print some statistics
        print "  %i functions found (%i parameters)\n" % (len(self.listOfFunctions), self.numberOfParams)

    def get_func_by_uname(self, uname):
        for func in self.listOfFunctions:
            if uname == func.uname:
                return func

    def create_input_file(self, fileName=None, fixAll=False):
        if fileName is None:
            fileName = "%s/results/temp_%s.dat" % (getcwd(), uuid.uuid4())
        fout = open(fileName, "w")
        fout.truncate(0)
        for func in self.listOfFunctions:
            fout.write(func.get_par_by_name("X0").tostring(fixAll))
            fout.write(func.get_par_by_name("Y0").tostring(fixAll))
            fout.write("FUNCTION " + func.name+"\n")
            for par in func.params[2:]:
                fout.write(par.tostring(fixAll))
        fout.close()
        # print "Model was saved to '%s'\n" % (fileName)
        return fileName

    def create_genome(self):
        """ This function creates a genome class for each
        parameter of the model. These classes will be stored in
        the genomeClasses list and their instances can be
        created using elements of this list (i.e. there is no
        names of these classes)"""
        genomeClasses = []
        genome = {}
        for func in self.listOfFunctions:
            for par in func.params:
                if not par.fixed:
                    class ParamGene(FloatGeneMax):
                        randMin = par.lowerLim
                        randMax = par.upperLim
                        mutProb = 0.5
                        mutAmt = 0.01
                    genomeClasses.append(ParamGene)
                    genome[func.uname+":"+par.name] = ParamGene
        return genome

    def genome_to_model(self, genome):
        """Set model parameter values according to given genome"""
        for gene in genome:
            funcName, parName = gene.split(":")
            func = self.get_func_by_uname(funcName)
            for par in func.params:
                if par.name == parName:
                    par.change_value(genome[gene])

    def check_boundaries(self, resModel):
        """ Method takes other model object and checks if its parameter
        values are close to the parameter limiths of the self model"""
        badParams = []
        for selfFunc, resFunc in zip(self.listOfFunctions, resModel.listOfFunctions):
            for selfParam, resParam in zip(selfFunc.params, resFunc.params):
                if selfParam.fixed:
                    # Fixed parameter, so it does not have boundaries. Nothing to check.
                    continue
                parRange = selfParam.upperLim - selfParam.lowerLim
                eps = parRange / 1000.0
                if (abs(resParam.value-selfParam.upperLim)<eps) or (abs(resParam.value-selfParam.lowerLim)<eps):
                    badParams.append("%s: %s" % (selfFunc.name, selfParam.name))
        return badParams

    def model_to_text(self, genNumber, textFile):
        """ Method saves current values of model parameters to a text file """
        if not exists(textFile):
            fout = open(textFile, "w", buffering=0)
            # Create a header as a first line of a file
            fout.write("# genNumber ")
            for func in self.listOfFunctions:
                for param in func.params:
                    fout.write("  %s.%s" % (func.uname, param.name))
            fout.write("\n")
        else:
            fout = open(textFile, "a", buffering=0)
        fout.write("%i " % genNumber)
        for func in self.listOfFunctions:
            for param in func.params:
                fout.write("  %9.3f" % param.value)
        fout.write("\n")


class GeneralParams(object):
    """ Input parameters unrelated to imfit (input image, size of the
    population etc.) """
    def __init__(self, fileName):
        for line in open(fileName):
            sLine = line.strip()
            if sLine.startswith("#"):
                # It is a comment line, just skip it
                continue
            if "#" in sLine:
                # Drop the comment part of the line if exists
                sLine = sLine[:sLine.index("#")].strip()
            if sLine.startswith("numOfCores"):
                self.numOfCores = int(sLine.split()[1])
                continue
            if sLine.startswith("Fits"):
                self.fitsToFit = sLine.split()[1]
                continue
            if sLine.startswith("PSF"):
                self.PSF = sLine.split()[1]
                continue
            if sLine.startswith("Mask"):
                self.mask = sLine.split()[1]
                continue
            if sLine.startswith("Weight"):
                self.weight = sLine.split()[1]
                continue
            if sLine.startswith("Popsize"):
                self.popSize = int(sLine.split()[1])
                continue
            if sLine.startswith("ZeroGenSize"):
                self.zeroGenSize = int(sLine.split()[1])
                continue
            if sLine.startswith("selectNbest"):
                self.selectNbest = int(sLine.split()[1])
                continue
            if sLine.startswith("addNew"):
                self.addNew = int(sLine.split()[1])
                continue
            if sLine.startswith("readNoise"):
                if sLine.split()[1] == "none":
                    self.readNoise == "none"
                else:
                    self.readNoise = float(sLine.split()[1])
                continue
            if sLine.startswith("gain"):
                if sLine.split()[1] == "none":
                    self.gain = "none"
                else:
                    self.gain = float(sLine.split()[1])
                continue
            if sLine.startswith("maxGenNumber"):
                self.maxGenNumber = int(sLine.split()[1])
                continue
            if sLine.startswith("fTol"):
                self.fTol = float(sLine.split()[1])
                continue
            if sLine.startswith("fSpan"):
                self.fSpan = int(sLine.split()[1])
                continue
            if sLine.startswith("imfitPath"):
                self.imfitPath = sLine.split()[1]
                continue
            if sLine.startswith("addImfitStr"):
                if sLine.split()[1] != "none":
                    self.addImfitStr = sLine.split('"')[1]
                else:
                    self.addImfitStr = " "
                continue
            if sLine.startswith("runLM"):
                self.runLM = sLine.split()[1]
                continue
            if sLine.startswith("numOfLM"):
                self.numOfLM = int(sLine.split()[1])
                continue
            if sLine.startswith("LMCores"):
                self.LMCores = int(sLine.split()[1])
                continue
            if sLine.startswith("saveGens"):
                self.saveGens = sLine.split()[1]
                continue
            if sLine.startswith("genTextFile"):
                if sLine.split()[1] == "none":
                    self.genTextFile == None
                else:
                    self.genTextFile = sLine.split()[1]
