# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
"""Tests for the geographic cache-key logic of the online weather providers.

These tests exercise ``_get_cache_filename`` in isolation, without performing
any network request, by building bare provider instances with
``object.__new__``. They guard against coordinate binning bugs where two
distinct locations on opposite sides of the equator or prime meridian would
collide onto the same cache key.
"""
import os
import unittest

from ..input.nasapower import NASAPowerWeatherDataProvider
from ..input.openmeteo import OpenMeteoWeatherDataProvider


class _CacheKeyTestMixin:
    """Shared assertions for provider cache-key binning.

    Subclasses must set ``provider_cls`` and implement ``make_filename``.
    """
    provider_cls = None

    def make_filename(self, lat, lon):
        """Return the cache filename for (lat, lon) without any network access."""
        # Build a bare instance so __init__ (which triggers a network request)
        # is never called; _get_cache_filename only reads a few attributes.
        obj = object.__new__(self.provider_cls)
        obj.model = "best_match"  # only used by the OpenMeteo provider
        return os.path.basename(obj._get_cache_filename(lat, lon))

    def test_no_cross_equator_collision(self):
        """Points just north and south of the equator must not share a cache key."""
        north = self.make_filename(0.05, 10.0)
        south = self.make_filename(-0.05, 10.0)
        self.assertNotEqual(north, south)

    def test_no_cross_meridian_collision(self):
        """Points just east and west of the prime meridian must not collide."""
        east = self.make_filename(52.0, 0.07)
        west = self.make_filename(52.0, -0.07)
        self.assertNotEqual(east, west)

    def test_uniform_binning_width(self):
        """The grid cell straddling zero must be 0.1 deg wide, like every other cell.

        With truncation toward zero the interval (-0.1, 0.1) collapses onto a
        single bin (double width). With flooring, -0.05 and 0.05 fall in
        adjacent, distinct bins.
        """
        below = self.make_filename(-0.05, 10.0)
        above = self.make_filename(0.05, 10.0)
        self.assertNotEqual(below, above)

    def test_same_cell_shares_key(self):
        """Two points inside the same 0.1 deg cell must share a cache key."""
        a = self.make_filename(52.51, 5.61)
        b = self.make_filename(52.59, 5.69)
        self.assertEqual(a, b)

    def test_northern_hemisphere_key_unchanged(self):
        """Positive coordinates keep their historical (int == floor) cache key."""
        fname = self.make_filename(52.56, 5.65)
        self.assertIn("LAT00525", fname)
        self.assertIn("LON00056", fname)


class Test_NASAPowerCacheKey(_CacheKeyTestMixin, unittest.TestCase):
    provider_cls = NASAPowerWeatherDataProvider


class Test_OpenMeteoCacheKey(_CacheKeyTestMixin, unittest.TestCase):
    provider_cls = OpenMeteoWeatherDataProvider


def suite():
    """This defines all the tests of a module"""
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Test_NASAPowerCacheKey))
    suite.addTest(unittest.makeSuite(Test_OpenMeteoCacheKey))
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
