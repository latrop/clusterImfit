#! /usr/bin/env python

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


class ImfitFunction(object):
    """ Class represents imfit function with
    all its parameters their ranges"""
    def __init__(self, funcName):
        self.funcName = funcName
        self.params = {}
    def add_parameter(self, newParameter):
        self.params[newParameter.name] = newParameter


class ImfitModel(object):
    """Imfit functions and their parameters"""
    def __init__(self, modelFileName):
        # Read imfit input file
        self.listOfFunctions = []
        funcName = None
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
                # If we are working already with some function, then
                # the list of parameters for this function is over and we can
                # add it to the function list
                if funcName is not None:
                    self.listOfFunctions.append(currentFunction)
                funcName = line.split()[1]
                currentFunction = ImfitFunction(funcName)
                currentFunction.add_parameter(x0)
                currentFunction.add_parameter(y0)
            else:
                # If line does not contain nor coordinates nor function name
                # then in has to be a parameter line
                param = parse_imfit_line(line)
                currentFunction.add_parameter(param)
