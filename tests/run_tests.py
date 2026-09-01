#!/usr/bin/env python3

"""Regression tests for GraviDy.

Every case runs a short simulation with one of the binaries and compares its
output against a reference stored in tests/reference.

The point of these tests is not to validate the physics, which the references
cannot do since they were produced by the code itself, but to notice when a
change to the code alters the result of a simulation.  Run them before and
after every refactoring.

The runs are done with a single OpenMP thread.  The energy of the system is
accumulated with «omp atomic» and the PN forces are written by more than one
thread, so with several threads the results are not reproducible.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# Colour escapes and timestamps of the log messages
ANSI_COLOUR = re.compile(r"\033\[[0-9;]*m")
TIMESTAMP = re.compile(r"\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\]")

# A number, in any of the notations that the code writes
NUMBER = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")

# How many differences of a case to print before giving up
MAX_REPORT = 10


def normalize(text, directories):
    """Remove everything that legitimately changes between two runs.

    That is, the colour escapes and the timestamps of the log messages, the
    ElapsedTime and the GFLOPS columns of the energy log, which are
    wall-clock measurements, and the directories of the run, which appear in
    the paths written to the .info file.
    """
    lines = []

    for line in text.splitlines():
        line = ANSI_COLOUR.sub("", line)
        line = TIMESTAMP.sub("[TIMESTAMP]", line)

        for directory, name in directories.items():
            line = line.replace(directory, name)

        fields = line.split()
        if len(fields) == 9 and (fields[0].isdigit() or
                                 (fields[0] == "#" and fields[1] == "IteTime")):
            line = " ".join(fields[:7])

        lines.append(line)

    return lines


def close_enough(expected, obtained, rtol, atol):
    """Compare two fields, with a tolerance when both of them are numbers.

    The tolerance is needed because the last digits of a floating point
    result depend on the compiler, on its optimisation flags and on the
    machine, while a real regression changes the numbers by much more.
    """
    if expected == obtained:
        return True

    if not (NUMBER.match(expected) and NUMBER.match(obtained)):
        return False

    a, b = float(expected), float(obtained)
    difference = abs(a - b)

    if difference <= atol:
        return True

    magnitude = max(abs(a), abs(b))

    return magnitude > 0 and difference / magnitude <= rtol


def differences(obtained, expected, rtol, atol):
    """List the differences between two normalized outputs."""
    found = []

    for number, (got, want) in enumerate(zip(obtained, expected), start=1):
        got_fields, want_fields = got.split(), want.split()

        if len(got_fields) != len(want_fields):
            found.append("line %d: %d fields instead of %d"
                         % (number, len(got_fields), len(want_fields)))
            continue

        for column, (a, b) in enumerate(zip(got_fields, want_fields), start=1):
            if not close_enough(b, a, rtol, atol):
                found.append("line %d: field %d is «%s» and not «%s»"
                             % (number, column, a, b))

    if len(obtained) > len(expected):
        found.append("%d lines more than the reference, which has %d"
                     % (len(obtained) - len(expected), len(expected)))
    elif len(obtained) < len(expected):
        found.append("%d lines are missing, the reference has %d"
                     % (len(expected) - len(obtained), len(expected)))

    return found


class Suite:
    """The set of cases, and the counters of what happened to them."""

    def __init__(self, options, work):
        self.bin_dir = options.bindir
        self.data_dir = Path(__file__).resolve().parent / "data"
        self.reference_dir = Path(__file__).resolve().parent / "reference"
        self.work = work
        self.rtol = options.rtol
        self.atol = options.atol
        self.update = options.update

        self.total = 0
        self.passed = 0
        self.failed = 0
        self.skipped = 0

        # References already written by --update in this run.  One shared by
        # two cases, like the one of the CPU and the MPI versions, is written
        # by the first of them and verified by the rest, otherwise the second
        # run would hide a difference between the two.
        self.written = set()

    def run(self, name, reference, program, arguments,
            launcher=(), outputs=None):
        """Run one case and compare its output against its reference.

        «outputs» are the files to compare, given as patterns relative to the
        working directory; without them the output of the program is used.
        """
        self.total += 1

        missing = self.missing_program(program, launcher)
        if missing is not None:
            print("  SKIP  %-16s (%s is not available)" % (name, missing))
            self.skipped += 1
            return

        log = self.work / (name + ".log")
        command = [*launcher, str(program), *arguments]

        environment = dict(os.environ)
        # A single thread, for the results to be reproducible
        environment["OMP_NUM_THREADS"] = "1"

        completed = subprocess.run(command, env=environment, cwd=self.work,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        log.write_bytes(completed.stdout)

        if completed.returncode != 0:
            print("  FAIL  %-16s (the run failed with status %d)"
                  % (name, completed.returncode))
            for line in completed.stdout.decode(errors="replace").splitlines():
                print("    " + line)
            self.failed += 1
            return

        if outputs is None:
            text = log.read_text(errors="replace")
        else:
            text = self.collect(outputs)

        self.compare(name, text, reference)

    def missing_program(self, program, launcher):
        """The first program of the case that cannot be run, if any."""
        if not (program.is_file() and os.access(program, os.X_OK)):
            return program.name

        if launcher and shutil.which(launcher[0]) is None:
            return launcher[0]

        return None

    def collect(self, patterns):
        """Join several output files into a single text."""
        text = ""

        for pattern in patterns:
            for path in sorted(self.work.glob(pattern)):
                text += "== %s ==\n" % path.name
                text += path.read_text(errors="replace")

        return text

    def compare(self, name, text, reference):
        directories = {str(self.work): "@WORK@", str(self.data_dir): "@DATA@"}
        obtained = normalize(text, directories)
        path = self.reference_dir / reference

        if self.update and reference not in self.written:
            path.write_text("".join(line + "\n" for line in obtained))
            self.written.add(reference)
            print("  WROTE %-16s (%s)" % (name, reference))
            self.passed += 1
            return

        if not path.is_file():
            print("  FAIL  %-16s (no reference, run with --update)" % name)
            self.failed += 1
            return

        expected = path.read_text().splitlines()
        found = differences(obtained, expected, self.rtol, self.atol)

        if not found:
            print("  ok    %-16s" % name)
            self.passed += 1
            return

        for message in found[:MAX_REPORT]:
            print("    " + message)
        if len(found) > MAX_REPORT:
            print("    ...")

        print("  FAIL  %-16s (differs from %s)" % (name, reference))
        self.failed += 1


def cases(suite):
    """Every case of the suite, in the order in which they are run."""
    cpu = suite.bin_dir / "gravidy-cpu"
    pn = suite.bin_dir / "gravidy-pn"
    mpi = suite.bin_dir / "gravidy-mpi"

    p128 = str(suite.data_dir / "p128.txt")
    two_body = str(suite.data_dir / "two-body.txt")

    # The energy log of a whole integration time
    suite.run("cpu-energy", "cpu-energy.ref", cpu,
              ["-i", p128, "-o", str(suite.work / "cpu-energy"),
               "-t", "1", "-p"])

    # The snapshots and the information file, which are the ones a user keeps
    suite.run("cpu-snapshots", "cpu-snapshots.ref", cpu,
              ["-i", p128, "-o", str(suite.work / "cpu-snapshots"),
               "-t", "0.25"],
              outputs=["cpu-snapshots.out.info",
                       "cpu-snapshots.out.snapshot_*"])

    # The Lagrange radii, which exercise the center of density and the sorting
    suite.run("cpu-lagrange", "cpu-lagrange.ref", cpu,
              ["-i", p128, "-o", str(suite.work / "cpu-lagrange"),
               "-t", "0.25", "-p", "-l"])

    # The state of every particle, the most detailed output of the code
    suite.run("cpu-all", "cpu-all.ref", cpu,
              ["-i", p128, "-o", str(suite.work / "cpu-all"),
               "-t", "0.25", "-p", "-a"])

    # The Post-Newtonian terms, which only accept two bodies
    suite.run("pn-two-body", "pn-two-body.ref", pn,
              ["-i", two_body, "-o", str(suite.work / "pn-two-body"),
               "-t", "0.5", "-p"])

    # The MPI version has to give the same numbers as the CPU one, so it is
    # compared against the very same reference
    suite.run("mpi-energy", "cpu-energy.ref", mpi,
              ["-i", p128, "-o", str(suite.work / "mpi-energy"),
               "-t", "1", "-p"],
              launcher=["mpirun", "-np", "2"])


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run the regression tests of GraviDy.",
        epilog="The cases whose binary is not built are skipped.")

    parser.add_argument("--update", action="store_true",
                        help="write the references instead of comparing "
                             "against them")
    parser.add_argument("--bindir", type=Path,
                        default=Path(__file__).resolve().parent.parent / "src",
                        help="where the binaries are (default: src)")
    parser.add_argument("--rtol", type=float, default=1e-9,
                        help="relative tolerance (default: 1e-9)")
    parser.add_argument("--atol", type=float, default=1e-14,
                        help="absolute tolerance (default: 1e-14)")

    options = parser.parse_args()

    if not options.bindir.is_dir():
        parser.error("%s does not exist" % options.bindir)
    options.bindir = options.bindir.resolve()

    return options


def main():
    options = parse_arguments()

    print("GraviDy regression tests")
    print("  binaries  %s" % options.bindir)
    print("  tolerance %g relative, %g absolute" % (options.rtol, options.atol))
    print("")

    with tempfile.TemporaryDirectory() as directory:
        suite = Suite(options, Path(directory))
        cases(suite)

    print("")
    print("%d cases: %d passed, %d failed, %d skipped"
          % (suite.total, suite.passed, suite.failed, suite.skipped))

    return 1 if suite.failed else 0


if __name__ == "__main__":
    sys.exit(main())
