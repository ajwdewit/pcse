# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014

import unittest
from datetime import date
import yaml
from ..agromanager import AgroManager
from ..base_classes import VariableKiosk
from ..exceptions import PCSEError


# Template for agromanagement testing
class TestAgroManagerSimpleTestTemplate(unittest.TestCase):

    agmt_input = None
    start_date = None
    end_date = None

    def runTest(self):
        d = yaml.load(self.agmt_input)
        kiosk = VariableKiosk()
        amgt = AgroManager(kiosk, d['AgroManagement'])
        self.assertEqual(amgt.start_date, self.start_date)
        self.assertEqual(amgt.end_date, self.end_date)

class TestAgroManager1(TestAgroManagerSimpleTestTemplate):

    agmt_input = """
                Version: 1.0
                AgroManagement:
                - 1999-08-01:
                    CropCalendar:
                        crop_id: winter-wheat
                        crop_start_date: 1999-09-15
                        crop_start_type: sowing
                        crop_end_date:
                        crop_end_type: maturity
                        max_duration: 300
                    TimedEvents:
                    StateEvents:
                - 2001-01-01:
                """
    start_date = date(1999, 8, 1)
    end_date = date(2001, 1, 1)


class TestAgroManager2(TestAgroManagerSimpleTestTemplate):

    agmt_input = """
                Version: 1.0
                AgroManagement:
                - 1999-08-01:
                    CropCalendar:
                        crop_id: winter-wheat
                        crop_start_date: 1999-09-15
                        crop_start_type: sowing
                        crop_end_date:
                        crop_end_type: maturity
                        max_duration: 300
                    TimedEvents:
                    StateEvents:
                """
    start_date = date(1999, 8, 1)
    end_date = date(2000, 7, 11)


class TestAgroManager3(TestAgroManagerSimpleTestTemplate):

    agmt_input = """
                Version: 1.0
                AgroManagement:
                - 1999-08-01:
                    CropCalendar:
                        crop_id: winter-wheat
                        crop_start_date: 1999-09-15
                        crop_start_type: sowing
                        crop_end_date: 2000-07-31
                        crop_end_type: harvest
                        max_duration: 300
                    TimedEvents:
                    StateEvents:
                """
    start_date = date(1999, 8, 1)
    end_date = date(2000, 7, 31)


class TestAgroManager4(TestAgroManagerSimpleTestTemplate):

    agmt_input = """
                Version: 1.0
                AgroManagement:
                - 1999-08-01:
                    CropCalendar:
                        crop_id: winter-wheat
                        crop_start_date: 1999-09-15
                        crop_start_type: sowing
                        crop_end_date: 2000-07-31
                        crop_end_type: harvest
                        max_duration: 300
                    TimedEvents:
                    -   event_signal: irrigate
                        name:  Timed irrigation events
                        comment: All irrigation amounts in cm
                        events_table:
                        - 2000-01-01: {irrigation_amount: 2, efficiency: 0.7}
                        - 2000-01-21: {irrigation_amount: 5, efficiency: 0.7}
                        - 2000-03-18: {irrigation_amount: 3, efficiency: 0.7}
                        - 2000-08-19: {irrigation_amount: 2, efficiency: 0.7}
                    StateEvents:
                """
    start_date = date(1999, 8, 1)
    end_date = date(2000, 8, 19)


class TestAgroManager5(TestAgroManagerSimpleTestTemplate):

    agmt_input = """
                Version: 1.0
                AgroManagement:
                - 1999-08-01:
                    CropCalendar:
                        crop_id: winter-wheat
                        crop_start_date: 1999-09-15
                        crop_start_type: sowing
                        crop_end_date:
                        crop_end_type: maturity
                        max_duration: 300
                    TimedEvents:
                    -   event_signal: irrigate
                        name:  Timed irrigation events
                        comment: All irrigation amounts in cm
                        events_table:
                        - 2000-01-01: {irrigation_amount: 2, efficiency: 0.7}
                        - 2000-01-21: {irrigation_amount: 5, efficiency: 0.7}
                        - 2000-03-18: {irrigation_amount: 3, efficiency: 0.7}
                        - 2000-03-19: {irrigation_amount: 2, efficiency: 0.7}
                    -   event_signal: apply_npk
                        name:  Timed N/P/K application table
                        comment: All fertilizer amounts in kg/ha
                        events_table:
                        - 2000-01-10: {N_amount : 10, P_amount: 5, K_amount: 2}
                        - 2000-01-31: {N_amount : 30, P_amount: 15, K_amount: 12}
                        - 2000-03-25: {N_amount : 50, P_amount: 25, K_amount: 22}
                        - 2000-04-05: {N_amount : 70, P_amount: 35, K_amount: 32}
                    StateEvents:
                - 2001-01-01:
                    CropCalendar:
                        crop_id: fodder-maize
                        crop_start_date: 2001-04-15
                        crop_start_type: sowing
                        crop_end_date:
                        crop_end_type: maturity
                        max_duration: 200
                    TimedEvents:
                    StateEvents:
                    -   event_signal: apply_npk
                        event_state: DVS
                        zero_condition: rising
                        name: DVS-based N/P/K application table
                        comment: all fertilizer amounts in kg/ha
                        events_table:
                        - 0.3: {N_amount : 1, P_amount: 3, K_amount: 4}
                        - 0.6: {N_amount: 11, P_amount: 13, K_amount: 14}
                        - 1.12: {N_amount: 21, P_amount: 23, K_amount: 24}
                    -   event_signal: irrigate
                        event_state: SM
                        zero_condition: falling
                        name: Soil moisture driven irrigation scheduling
                        comment: all irrigation amounts in cm
                        events_table:
                        - 0.15: {irrigation_amount: 2, efficiency: 0.7}
                - 2001-12-15:
                """

    start_date = date(1999, 8, 1)
    end_date = date(2001, 12, 15)


# This should raise an error because the last campaign has StateEvents but no explicit end date
# is given.
class TestAgroManager6(unittest.TestCase):
    agmt_input = """
                Version: 1.0
                AgroManagement:
                - 2001-01-01:
                    CropCalendar:
                    TimedEvents:
                    StateEvents:
                    -   event_signal: irrigate
                        event_state: SM
                        zero_condition: falling
                        name: Soil moisture driven irrigation scheduling
                        comment: all irrigation amounts in cm
                        events_table:
                        - 0.15: {irrigation_amount: 2, efficiency: 0.7}
                """

    def runTest(self):
        d = yaml.load(self.agmt_input)
        kiosk = VariableKiosk()
        self.amgt = AgroManager(kiosk, d['AgroManagement'])
        self.assertRaises(PCSEError, self._get_end_date)

    def _get_end_date(self):
        return self.amgt.end_date


# This should raise an error because the timed event on 2001-03-18 is not within the
# campaign interval
class TestAgroManager7(unittest.TestCase):
    agmt_input = """
                Version: 1.0
                AgroManagement:
                - 2000-01-01:
                    CropCalendar:
                    TimedEvents:
                     - event_signal: irrigate
                       name:  Timed irrigation events
                       comment: All irrigation amounts in cm
                       events_table:
                       - 2000-01-01: {irrigation_amount: 2, efficiency: 0.7}
                       - 2000-01-21: {irrigation_amount: 5, efficiency: 0.7}
                       - 2001-03-18: {irrigation_amount: 3, efficiency: 0.7}
                       - 2000-03-19: {irrigation_amount: 2, efficiency: 0.7}
                    StateEvents:
                - 2001-01-01: null
                """

    def runTest(self):
        self.assertRaises(PCSEError, self._start_agromanager)

    def _start_agromanager(self):
        d = yaml.load(self.agmt_input)
        kiosk = VariableKiosk()
        amgt = AgroManager(kiosk, d['AgroManagement'])


# This should raise an error because the crop start date is not within the
# campaign interval
class TestAgroManager8(unittest.TestCase):
    agmt_input = """
                Version: 1.0
                AgroManagement:
                - 2000-01-01:
                    CropCalendar:
                        crop_id: fodder-maize
                        crop_start_date: 2001-04-15
                        crop_start_type: sowing
                        crop_end_date:
                        crop_end_type: maturity
                        max_duration: 200
                    TimedEvents:
                    StateEvents:
                - 2001-01-01: null
                """

    def runTest(self):
        self.assertRaises(PCSEError, self._start_agromanager)

    def _start_agromanager(self):
        d = yaml.load(self.agmt_input)
        kiosk = VariableKiosk()
        amgt = AgroManager(kiosk, d['AgroManagement'])

def suite():
    """ This defines all the tests of a module"""
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestAgroManager1))
    suite.addTest(unittest.makeSuite(TestAgroManager2))
    suite.addTest(unittest.makeSuite(TestAgroManager3))
    suite.addTest(unittest.makeSuite(TestAgroManager4))
    suite.addTest(unittest.makeSuite(TestAgroManager5))
    suite.addTest(unittest.makeSuite(TestAgroManager6))
    suite.addTest(unittest.makeSuite(TestAgroManager7))
    suite.addTest(unittest.makeSuite(TestAgroManager8))
    return suite

if __name__ == '__main__':
   unittest.TextTestRunner(verbosity=2).run(suite())