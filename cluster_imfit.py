#! /usr/bin/env python

from libs.read_input import ImfitModel


model = ImfitModel("input.dat")

model.create_input_file("out.dat")
