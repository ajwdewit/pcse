# -*- coding: utf-8 -*-
# Copyright (c) 2004-2018 Wageningen Environmental Sciences, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), January 2018
"""This allows to run the YAML test files using `python -m tests`, optionally enabling the
full test suite with `python -m tests --full`
"""
import argparse
import unittest

from .run_tests import make_test_suite


def create_parser():
    parser = argparse.ArgumentParser(description='Run PCSE test suite', prog="python -m tests")
    parser.add_argument('--full', dest='full_tests', action='store_true', default=False,
                        help='Run full test suite instead of short suite', required=False,
                        )
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    if args.full_tests is True:
        suite = make_test_suite(quick=False)
    else:
        suite = make_test_suite(quick=True)

    unittest.TextTestRunner(verbosity=2).run(suite)


if __name__ == "__main__":
    main()