#!/usr/bin/env bash
# Cell Quant Studio Universal - run Jython scripts against a Fiji installation
# from a terminal (Linux, macOS, or Git Bash on Windows).
#
#   tests/run_tests.sh /path/to/Fiji                 runs tests/run_tests.py
#   tests/run_tests.sh /path/to/Fiji some_script.py  runs the given script(s)
#
# The Fiji folder can also be given in the FIJI_DIR environment variable.
# The script uses the Java runtime bundled with Fiji when there is one, and
# the system "java" otherwise. Exit status 0 means every test passed.

set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
FIJI="${FIJI_DIR:-}"
if [ $# -gt 0 ] && [ -d "$1/jars" ]; then
    FIJI="$1"
    shift
fi
if [ -z "$FIJI" ] || [ ! -d "$FIJI/jars" ]; then
    echo "Give the Fiji folder (the one holding 'jars' and 'plugins') as the first argument or in FIJI_DIR." >&2
    exit 2
fi

JAVA="$(find "$FIJI/java" -type f \( -name java -o -name java.exe \) 2>/dev/null | head -n 1)"
if [ -z "$JAVA" ]; then
    JAVA="java"
fi

SEP=":"
case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*) SEP=";" ;;
esac

if [ $# -eq 0 ]; then
    set -- "$HERE/run_tests.py"
fi

cd "$FIJI" || exit 2
# Extra Java options, for example CELLQUANT_JAVA_OPTS="-Xmx8g" for large images.
# shellcheck disable=SC2086
exec "$JAVA" ${CELLQUANT_JAVA_OPTS:-} -Dpython.cachedir.skip=true -Dpython.console.encoding=UTF-8 \
    -Dplugins.dir="$FIJI" \
    -cp "jars/*${SEP}jars/bio-formats/*${SEP}plugins/*" \
    org.python.util.jython "$HERE/cli_runner.py" "$@"
