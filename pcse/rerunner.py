from __future__ import print_function
import os
from .exceptions import PCSEError
import textwrap

class Rerunner(list):
    """Class for reading and iterating through an FSE rerun file.

    Rerun files are used for executing many model runs with varying inputs.
    This was particularly useful for FSE models that were not easily
    scriptable and adjustable. For PCSE models this functionality is
    less relevant as the model can be easily embedded within a python
    loop to iterate over combinations of input values.

    Nevertheless, as the concept of rerun files is still highly popular
    in Wageningen this class was added to accomodate users that like to
    work with rerun files.

    This class parses a rerun file and stores the individual reruns
    as a list of dictionaries. Each new rerun in the file should be
    separated by a line stating "RUNNAME='<run name>';"

     :

        >>>reruns = Rerunner('reruns.txt')
        >>>for rerun in reruns:
            # processing goes here..
    """

    def __init__(self, rerun_file):

        try:
            self.rerun_file_fp = os.path.abspath(rerun_file)
            with open(self.rerun_file_fp) as fp:
                lines = fp.readlines()
        except IOError as exc:
            msg = "Failed opening '%s': %s" % (self.rerun_file_fp, exc.message)
            raise PCSEError(msg)

        rerun_set = None
        try:
            for i, line in enumerate(lines):
                ln = i+1
                line = line.strip() #.replace(" ","")
                if not line:  # empty line
                    continue
                if line.startswith("*"):
                    continue
                # if not line.endswith(";"):
                #     msg = "';' character missing on line %i" % ln
                #     raise PCSEError(msg)

                name, strvalue = line.split("=", 1)
                name = name.strip()
                value = eval(strvalue)

                if name.upper() == "RUNNAM":
                    if not rerun_set is None:
                        self.append(rerun_set)
                    rerun_set = {name: value}
                else:
                    rerun_set[name] = value

            # Append last rerun set
            if rerun_set is not None:
                self.append(rerun_set)

        except (SyntaxError, ValueError) as exc:
            msg = "Failed parsing line %i." % ln
            raise PCSEError(msg)

    def __str__(self):

        msg = "Rerun loaded from: %s\n" % self.rerun_file_fp
        if len(self) > 0:
            msg += ("Containing %i rerun sets\n" % len(self))
            d = "\n".join(textwrap.wrap("First rerun set: %s\n" % self[0], subsequent_indent="    "))
            msg += (d + "\n")
            d = "\n".join(textwrap.wrap("Last rerun set: %s\n" % self[-1], subsequent_indent="    "))
            msg += d
        else:
            msg += "File contains no rerun sets.\n"
        return msg
