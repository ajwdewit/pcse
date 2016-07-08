__author__ = 'wit015'
import sys, os
import yaml
from .. import exceptions as exc

class YAMLAgroManagementReader(list):
    """Reads PCSE agromanagement files in the YAML format.

    :param fname: filename of the agromanagement file. If fname is not provided as a absolute or
        relative path the file is assumed to be in the current working directory.
    """

    def __init__(self, fname):
        fname_fp = os.path.normpath(os.path.abspath(fname))
        if not os.path.exists(fname_fp):
            msg = "Cannot find agromanagement file: %s" % fname_fp
            raise exc.PCSEError(msg)

        with open(fname) as fp:
            try:
                r = yaml.load(fp)
            except yaml.YAMLError as e:
                msg = "Failed parsing agromanagement file %s: %s" % (fname_fp, e)
                raise exc.PCSEError(msg)

        list.__init__(self, r['AgroManagement'])

    def __str__(self):
        return yaml.dump(self, default_flow_style=False)
