#!/bin/sh
# Compiles TotalRepeats together with its test suite and runs the suite.
# Requires only a JDK (25+); there are no external test dependencies.
#
#   ./run-tests.sh
#
# Exits non-zero if any check fails, so it can gate a build or a CI job.
set -e

here=$(cd "$(dirname "$0")" && pwd)
out="$here/build/test-classes"

rm -rf "$out"
mkdir -p "$out"

javac -d "$out" "$here"/src/*.java "$here"/tests/*.java

exec java -cp "$out" TestRunner
