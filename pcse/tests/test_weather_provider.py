# -*- coding: utf-8 -*-
import unittest
from datetime import date, timedelta

from ..base.weather import WeatherDataProvider


class EnsembleWeatherDataProvider(WeatherDataProvider):
    supports_ensembles = True


class TestWeatherDataProvider(unittest.TestCase):
    def setUp(self):
        self.provider = EnsembleWeatherDataProvider()
        self.start = date(2000, 1, 1)

    def add_ensemble_day(self, day):
        for member_id in range(3):
            self.provider.store[(day, member_id)] = object()

    def test_ensemble_members_do_not_count_as_separate_days(self):
        for offset in range(3):
            self.add_ensemble_day(self.start + timedelta(days=offset))

        self.assertEqual(self.provider.missing, 0)

    def test_missing_count_uses_unique_days_when_ensemble_day_is_absent(self):
        self.add_ensemble_day(self.start)
        self.add_ensemble_day(self.start + timedelta(days=2))

        self.assertEqual(self.provider.missing, 1)


def suite():
    return unittest.defaultTestLoader.loadTestsFromTestCase(TestWeatherDataProvider)


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
