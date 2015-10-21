#! /usr/bin/env python


from pygene.gene import FloatGene, FloatGeneMax
from pygene.organism import Organism, MendelOrganism
from pygene.population import Population

from libs.read_input import ImfitModel


def create_Converger_class(model):
    class Converger(MendelOrganism):
        """
        Implements the organism which tries to converge a function
        """
        model = model
        genome = model.create_genome()

        def fitness(self):
            """
            Implements the 'fitness function' for this species.
            Organisms try to evolve to minimise this function's value
            """
            params = {}
            for key in self.__dict__['genes'].keys():
                params[key] = self[key]
            self.model.genome_to_model(params)
    return Converger


model = ImfitModel("input.dat")
Converger = create_Converger_class(model)
c = Converger()
c.fitness()
