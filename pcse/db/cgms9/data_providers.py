# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2015
"""
Data providers for weather, agromanagement, soil, crop and site data. Also
a class for testing STU suitability for a given crop.

Data providers are compatible with a CGMS 9 database schema.
"""
import datetime

from sqlalchemy import MetaData, select, Table, and_
import yaml

from ...util import check_date
from ... import exceptions as exc
from ..cgms11.data_providers import fetch_crop_name


class AgroManagementDataProvider(list):
    """Class for providing agromanagement data from the CROP_CALENDAR table in a CGMS9 database.

    :param engine: SqlAlchemy engine object providing DB access
    :param grid_no: Integer grid ID, maps to the grid_no column in the table
    :param crop_no: Integer id of crop, maps to the crop_no column in the table
    :param campaign_year: Integer campaign year, maps to the YEAR column in the table.
           The campaign year refers to the year of the crop start. Thus for crops
           crossing calendar years, the start_date can be in the previous year as the
           harvest.

    Note that by default the campaign_start_date is set equal to the crop_start_date which
    means that the simulation starts when the crop starts. In some cases this is undesirable
    and an earlier start date should be used. In that case use the `set_campaign_start_date(date)`
    to update the campaign_start_date.
    """
    agro_management_template = """
          - {campaign_start_date}:
                CropCalendar:
                    crop_id: '{crop_name}'
                    crop_start_date: {crop_start_date}
                    crop_start_type: {crop_start_type}
                    crop_end_date: {crop_end_date}
                    crop_end_type: {crop_end_type}
                    max_duration: {duration}
                TimedEvents: null
                StateEvents: null
          - {campaign_end_date}: null
        """

    def __init__(self, engine, grid_no, crop_no, campaign_year):
        list.__init__(self)
        self.grid_no = int(grid_no)
        self.crop_no = int(crop_no)
        self.campaign_year = int(campaign_year)
        self.crop_name = fetch_crop_name(engine, self.crop_no)
        self.db_resource = str(engine)[7:-1]

        metadata = MetaData(engine)
        table_cc = Table("crop_calendar", metadata, autoload=True)

        r = select([table_cc], and_(table_cc.c.grid_no == self.grid_no,
                                    table_cc.c.crop_no == self.crop_no,
                                    table_cc.c.year == self.campaign_year)).execute()
        row = r.fetchone()
        r.close()
        if row is None:
            msg = "Failed deriving crop calendar for grid_no %s, crop_no %s " % (grid_no, crop_no)
            raise exc.PCSEError(msg)

        # Determine the start date.
        year = int(row.year)
        month = int(row.start_month1)
        day = int(row.start_monthday1)
        self.crop_start_date = check_date(datetime.date(year, month, day))
        self.campaign_start_date = self.crop_start_date

        # Determine the start date/type. Only sowing|emergence is accepted by PCSE/WOFOST
        cgms_start_type = str(row.start_type).strip()
        if cgms_start_type == "fixed_sowing":
            self.crop_start_type = "sowing"
        elif cgms_start_type == "fixed_emergence":
            self.crop_start_type = "emergence"
        else:
            msg = "Unsupported START_TYPE in CROP_CALENDAR table: %s" % row.start_type
            raise exc.PCSEError(msg)

        # Determine maximum duration of the crop
        self.max_duration = int(row.max_duration)

        # Determine crop end date/type and the end of the campaign
        self.crop_end_type = str(row.end_type).strip().lower()
        if self.crop_end_type not in ["harvest", "earliest", "maturity"]:
            msg = ("Unrecognized option for END_TYPE in table "
                   "CROP_CALENDAR: %s" % row.end_type)
            raise exc.PCSEError(msg)

        if self.crop_end_type == "maturity":
            self.crop_end_date = "null"
            self.campaign_end_date = self.crop_start_date + datetime.timedelta(days=self.max_duration)
        else:
            month = int(row.end_month)
            day = int(row.end_monthday)
            self.crop_end_date = datetime.date(year, month, day)
            if self.crop_end_date <= self.crop_start_date:
                self.crop_end_date = datetime.date(year+1, month, day)
            self.campaign_end_date = self.crop_end_date

        input = self._build_yaml_agromanagement()
        self._parse_yaml(input)

    def _build_yaml_agromanagement(self):
        """Builds the YAML agromanagent string"""

        input = self.agro_management_template.format(campaign_start_date=self.campaign_start_date,
                                                     crop_name=self.crop_name,
                                                     crop_start_date=self.crop_start_date,
                                                     crop_start_type=self.crop_start_type,
                                                     crop_end_date=self.crop_end_date,
                                                     crop_end_type=self.crop_end_type,
                                                     duration=self.max_duration,
                                                     campaign_end_date=self.campaign_end_date
                                                     )
        return input

    def _parse_yaml(self, input):
        """Parses the input YAML string and assigns to self"""
        try:
            items = yaml.load(input)
        except yaml.YAMLError as e:
            msg = "Failed parsing agromanagement string %s: %s" % (input, e)
            raise exc.PCSEError(msg)
        self.extend(items)

    def set_campaign_start_date(self, start_date):
        """Updates the value for the campaign_start_date.

        This is useful only when the INITIAL_SOIL_WATER table in CGMS9 defines a different
        campaign start date which should be used instead.
        """
        self.campaign_start_date = check_date(start_date)
        input = self._build_yaml_agromanagement()
        self._parse_yaml(input)
