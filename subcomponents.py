#!/usr/bin/env python

import argparse
import subprocess
import os
from libs.read_input import ImfitModel


def main(args):
    model = ImfitModel(args.imfit)
    for func in model.listOfFunctions:
        # Create an imfit file
        imfit_model_name = "%s.imfit" % func.uname
        func.to_model_file(imfit_model_name)
        # Call imfit to create a fits file
        fits_name = "%s.fits" % func.uname
        runString = "makeimage %s --refimage %s --output %s" % (imfit_model_name, args.fits, fits_name)
        proc = subprocess.Popen(runString, shell=True)
        proc.wait()
        os.remove(imfit_model_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("imfit", type=str,
                        help="File with imfir model.")
    parser.add_argument("fits", type=str,
                        help="File with FITS file")
    args = parser.parse_args()
    print(args)
    main(args)
