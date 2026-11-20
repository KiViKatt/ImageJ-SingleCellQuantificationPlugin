# Cell Quant Studio Universal - command-line runner for Jython scripts
#
# Runs one or more Jython scripts with Fiji's classes on the class path and
# ends the Java process with an exit status: 0 when every script finished,
# 1 when one of them raised. It exists because a plain
#   java ... org.python.util.jython tests/run_tests.py
# does not end by itself once ImageJ has started its window toolkit.
#
# It is started by tests/run_tests.sh; see that file. (On Windows, run that
# script from Git Bash.)
#
# The source is kept plain ASCII with no encoding declaration on purpose.

import sys
import traceback

from java.lang import System
from java.lang import Throwable as JavaThrowable

# Fiji reads its saved options when it starts. Without a running Fiji they are
# read lazily, in the middle of the first file-open, which would look like a
# script changing a global option. Read them now instead.
try:
    from ij import Prefs
    Prefs.load(None, None)
except (Exception, JavaThrowable):
    pass

cli_status = 0
for cli_path in sys.argv[1:]:
    cli_namespace = {"__file__": cli_path, "__name__": "__main__"}
    try:
        execfile(cli_path, cli_namespace)
    except SystemExit as cli_exit:
        if cli_exit.code not in (None, 0):
            cli_status = 1
    except (Exception, JavaThrowable):
        traceback.print_exc()
        cli_status = 1

sys.stdout.flush()
sys.stderr.flush()
System.exit(cli_status)
