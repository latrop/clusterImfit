#! /usr/bin/env python

from pygene.gene import FloatGene, FloatGeneMax, IntGene
from pygene.organism import Organism, MendelOrganism


def parse_imfit_line(line):
    """ Function parses line of imfit data file and
    returns parameters"""
    params = line.split()
    name = params[0]
    value = float(params[1]) # Value of the parameter is the second entry of the line
    # (after the name)
    # Lets now find range of values. Its have to contain the coma and mast
    # be the second enrty (otherwise imfin wont work)
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
    def tostring(self):
        return "{:12}{:6.2f}{:12.2f},{:1.2f}\n".format(self.name, self.value, self.lowerLim, self.upperLim)
    def change_value(self, newValue):
        self.value = newValue
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
        funcName = None
        ident = -1
        for line in open(modelFileName):
            sLine = line.strip()
            if line.startswith("#"):
                # It is a comment line, just skip it
                continue
            if "#" in line:
                # Drop the comment part of the line if exists
                line = line[:line.index("#")].strip()
            if line.startswith("X0"):
                x0 = parse_imfit_line(line)
            elif line.startswith("Y0"):
                y0 = parse_imfit_line(line)
            elif line.startswith("FUNCTION"):
                # New function is found.
                ident += 1
                # If we are working already with some function, then
                # the list of parameters for this function is over and we can
                # add it to the function list
                if funcName is not None:
                    self.listOfFunctions.append(currentFunction)
                funcName = line.split()[1]
                currentFunction = ImfitFunction(funcName, ident)
                currentFunction.add_parameter(x0)
                currentFunction.add_parameter(y0)
            else:
                # If line does not contain nor coordinates nor function name
                # then in has to be a parameter line
                param = parse_imfit_line(line)
                currentFunction.add_parameter(param)
        # append the last function
        self.listOfFunctions.append(currentFunction)
        # Print some statistics
        print "  %i functions found\n" % len(self.listOfFunctions)

    def get_func_by_uname(self, uname):
        for func in self.listOfFunctions:
            if uname == func.uname:
                return func

    def create_input_file(self, fileName):
        fout = open(fileName, "w")
        fout.truncate(0)
        for func in self.listOfFunctions:
            fout.write(func.get_par_by_name("X0").tostring())
            fout.write(func.get_par_by_name("Y0").tostring())
            fout.write(func.name+"\n")
            for par in func.params[2:]:
                fout.write(par.tostring())
        fout.close()
        print "Model was saved to '%s'\n" % (fileName)

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
                class ParamGene(FloatGeneMax):
                    randMin = par.lowerLim
                    randMax = par.upperLim
                    mutProb = 0.5
                    mutAmt = 0.05
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
        self.create_input_file("run.dat")
