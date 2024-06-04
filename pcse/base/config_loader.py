import os
from pathlib import Path
from .. import exceptions as exc
import textwrap


class ConfigurationLoader(object):
    """Class for loading the model configuration from a PCSE configuration files

        :param config: string given file name containing model configuration
        """
    _required_attr = ("CROP", "SOIL", "AGROMANAGEMENT", "OUTPUT_VARS", "OUTPUT_INTERVAL",
                      "OUTPUT_INTERVAL_DAYS", "SUMMARY_OUTPUT_VARS")
    defined_attr = []
    model_config_file = None
    description = None

    def __init__(self, config):

        if not isinstance(config, (str, Path)):
            msg = ("Keyword 'config' should provide the name of the file (string or pathlib.Path)" +
                   "storing the configuration of the model PCSE should run.")
            raise exc.PCSEError(msg)

        # check if model configuration file is an absolute or relative path. If
        # not assume that it is located in the 'conf/' folder in the PCSE
        # distribution
        config = Path(config)
        if config.is_absolute():
            mconf = config
        else:
            this_dir = Path(__file__).parent
            pcse_dir = this_dir.parent
            mconf = pcse_dir / "conf" / config
        model_config_file = mconf.resolve()

        # check that configuration file exists
        if not model_config_file.exists():
            msg = "PCSE model configuration file does not exist: %s" % model_config_file
            raise exc.PCSEError(msg)
        # store for later use
        self.model_config_file = model_config_file

        # Load file using execfile
        try:
            loc = {}
            bytecode = compile(open(model_config_file).read(), model_config_file, 'exec')
            exec(bytecode, {}, loc)
        except Exception as e:
            msg = "Failed to load configuration from file '%s' due to: %s"
            msg = msg % (model_config_file, e)
            raise exc.PCSEError(msg)

        # Add the descriptive header for later use
        if "__doc__" in loc:
            desc = loc.pop("__doc__")
            if len(desc) > 0:
                self.description = desc
                if self.description[-1] != "\n":
                    self.description += "\n"

        # Loop through the attributes in the configuration file
        for key, value in list(loc.items()):
            if key.isupper():
                self.defined_attr.append(key)
                setattr(self, key, value)

        # Check for any missing compulsary attributes
        req = set(self._required_attr)
        diff = req.difference(set(self.defined_attr))
        if diff:
            msg = "One or more compulsory configuration items missing: %s" % list(diff)
            raise exc.PCSEError(msg)

    def __str__(self):
        msg = "PCSE ConfigurationLoader from file:\n"
        msg += "  %s\n\n" % self.model_config_file
        if self.description is not None:
            msg += ("%s Header of configuration file %s\n"% ("-"*20, "-"*20))
            msg += self.description
            if msg[-1] != "\n":
                msg += "\n"
            msg += ("%s Contents of configuration file %s\n"% ("-"*19, "-"*19))
        for k in self.defined_attr:
             r = "%s: %s" % (k, getattr(self, k))
             msg += (textwrap.fill(r, subsequent_indent="  ") + "\n")
        return msg

    def update_output_variable_lists(self, output_vars=None, summary_vars=None, terminal_vars=None):
        """Updates the lists of output variables that are defined in the configuration file.

        This is useful because sometimes you want the flexibility to get access to an additional model
        variable which is not in the standard list of variables defined in the model configuration file.
        The more elegant way is to define your own configuration file, but this adds some flexibility
        particularly for use in jupyter notebooks and exploratory analysis.

        Note that there is a different behaviour given the type of the variable provided. List and string
        inputs will extend the list of variables, while set/tuple inputs will replace the current list.

        :param output_vars: the variable names to add/replace for the OUTPUT_VARS configuration variable
        :param summary_vars: the variable names to add/replace for the SUMMARY_OUTPUT_VARS configuration variable
        :param terminal_vars: the variable names to add/replace for the TERMINAL_OUTPUT_VARS configuration variable
        """
        config_varnames = ["OUTPUT_VARS", "SUMMARY_OUTPUT_VARS", "TERMINAL_OUTPUT_VARS"]
        for varitems, config_varname in zip([output_vars, summary_vars, terminal_vars], config_varnames):
            if varitems is None:
                continue
            else:
                if isinstance(varitems, str):  # A string: we extend the current list
                    getattr(self, config_varname).extend(varitems.split())
                elif isinstance(varitems, list):  # a list: we extend the current list
                    getattr(self, config_varname).extend(varitems)
                elif isinstance(varitems, (tuple, set)):  # tuple/set we replace the current list
                    setattr(self, config_varname, list(varitems))
                else:
                    msg = "Unrecognized input for `output_vars` to engine(): %s" % output_vars
                    print(msg)
