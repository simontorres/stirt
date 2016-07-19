#!/usr/local/python
"""CCD Reduction Routine for e2v camera at Swope Telescope

This pipeline performs all the necessary reduction processes for single nights

"""

import os
from astropy.io import fits
import ccdproc as ccd
import pandas as pd
import argparse
import logging as log
import warnings

warnings.filterwarnings('ignore')
log.basicConfig(level=log.INFO)

__author__ = 'Simon Torres'
__date__ = '2016-06-28'
__version__ = "0.1"
__email__ = "storres@ctio.noao.edu"
__status__ = "Development"


class MainApp:
    """docstrings

    """

    def __init__(self):
        self.args = self.get_args()

    def get_args(self):
        """Handles the argparse library and returns the arguments



        Returns:

        """

        parser = argparse.ArgumentParser(description='CCD image reduction tool for the Enrietta Swope Telescope.')

        parser.add_argument('-p', '--data-path',
                            action='store',
                            default='./',
                            type=str,
                            metavar='<Source Path>',
                            dest='source',
                            help='Path for location of raw data. Default <./>')

        parser.add_argument('-d', '--proc-path',
                            action='store',
                            default='./',
                            type=str,
                            metavar='<Destination Path>',
                            dest='destiny',
                            help='Path for destination of processed data. Default <./>')

        parser.add_argument('-v', '--debug',
                            action='store_true',
                            default=False,
                            dest='debug',
                            help='Enable debug mode')

        parser.add_argument('-sc', '--shutter',
                            action='store_true',
                            default=False,
                            dest='shutter',
                            help='Apply shutter correction')

        parser.add_argument('--skip-overscan',
                            action='store_false',
                            default=True,
                            dest='overscan',
                            help='Skip overscan correction and trim process.\
                                 Helpful when you want to repeat later processes')

        parser.add_argument('--skip-bias',
                            action='store_false',
                            default=True,
                            dest='bias',
                            help='Skip bias correction. Helpful when you want to repeat later processes')

        parser.add_argument('--skip-linearity',
                            action='store_false',
                            default=True,
                            dest='linearity',
                            help='Skip linearity correction. Helpful when you want to repeat later processes')

        parser.add_argument('--skip-flats',
                            action='store_false',
                            default=True,
                            dest='flats',
                            help='Skip flat correction. Helpful when you want to repeat later processes')

        parser.add_argument('--single-chip',
                            action='store',
                            default=5,
                            type=int,
                            metavar='<Chip Number>',
                            dest='chipNumber',
                            help='Process single chip instead of all. Default all')

        parser.add_argument('--mosaic',
                            action='store_true',
                            default=False,
                            dest='mosaic',
                            help='Build mosaic')
        parser.add_argument('--clean',
                            action='store_true',
                            default=False,
                            dest='clean',
                            help='Clean directory')

        args = parser.parse_args()
        if not os.path.isdir(args.source):
            log.error("Source Directory Doesn't exist.")
            parser.exit("Bye!")
        if not os.path.isdir(args.destiny):
            log.warning("Destination folder doesn't exist")
            if os.getcwd() not in args.destiny:
                log.warning("Please Create the destination folder: %s"%args.destiny)
                parser.exit("Bye!")
            else:
                log.info("Creating Destination folder: %s"%args.destiny)
                os.makedirs(args.destiny)

        return args


if __name__ == '__main__':
    App = MainApp()