#!/usr/bin/env python

import os
import shutil
import sys
import getopt

from rtl2 import test_util
from rtl2 import params
from rtl2 import test_report as report


def reg_test_gc(argv):
    usage = """
    ./reg_test_gc [--before|-b 2000-00-00]
       testfile.ini
    """

    if len(sys.argv) == 1:
        print(usage)
        sys.exit(2)

    try:
        opts, next_ = getopt.getopt(argv[1:], "b:", ["before="])

    except getopt.GetoptError:
        print("invalid calling sequence")
        print(usage)
        sys.exit(2)

    # defaults
    gcdate = ""

    for o, a in opts:
        if o in ["--before", "-b"]:
            gcdate = a

    try:
        testFile = next_[0]

    except IndexError:
        print("ERROR: a test file was not specified")
        print(usage)
        sys.exit(2)

    if not gcdate:
        print("ERROR: date was not specified")
        print(usage)
        sys.exit(2)

    gcd = valid_date(gcdate)
    if gcd == "":
        print("ERROR: invalid date", gcdate)
        print(usage)
        sys.exit(2)

    print("loading ", testFile)

    args = test_util.get_args([testFile])

    suite, testList = params.load_params(args)
    activeTestList = [t.name for t in testList]

    benchmarkTestList = [t for t in testList if not (t.compileTest or t.restartTest)]
    benchmarkNotFound = {}
    for t in benchmarkTestList:
        benchmarkNotFound[t.name] = ""

    ### clean up the web dir
    print("\ncleaning ", suite.webTopDir)

    os.chdir(suite.webTopDir)
    validDirs = []
    for d in os.listdir(suite.webTopDir):
        if d.startswith("20") and os.path.isdir(d):
            statusFile = d + "/" + d + ".status"
            if os.path.isfile(statusFile):
                validDirs.append(d)
    validDirs.sort()
    validDirs.reverse()

    latestBMDate = {}

    for d in validDirs:
        bmtests = list(benchmarkNotFound)
        if d >= gcd and bmtests:
            if isBenchmarkDir(d):
                for t in bmtests:
                    if findBenchmark(d, t):
                        del benchmarkNotFound[t]
                        latestBMDate[t] = d
        else:
            if isBenchmarkDir(d) and bmtests:
                found = False
                for t in bmtests:
                    if findBenchmark(d, t):
                        found = True
                        del benchmarkNotFound[t]
                        latestBMDate[t] = d
                if not found:
                    rmDir(d)
            else:
                rmDir(d)

    ### clean up the test dir
    testDirs = os.path.join(suite.testTopDir, suite.suiteName + "-tests")
    print("\ncleaning ", testDirs)

    os.chdir(testDirs)
    validDirs = []
    for d in os.listdir(testDirs):
        if d.startswith("20") and os.path.isdir(d):
            validDirs.append(d)
    validDirs.sort()
    validDirs.reverse()

    for d in validDirs:
        if d < gcd:
            tests = [t for t in os.listdir(d) if os.path.isdir(os.path.join(d, t))]
            found = False
            for t in tests:
                if t in latestBMDate.keys() and latestBMDate[t] == d:
                    found = True
                    break
            if not found:
                rmDir(d)

    print("\ncreating suite report...")
    report.report_all_runs(suite, activeTestList)

    print("\nGarbage cleaning finished.")


def valid_date(gcdate):
    try:
        y, m, d = gcdate.split("-")
        yi = int(y)
        mi = int(m)
        di = int(d)
        return (
            "-".join([y, m.zfill(2), d.zfill(2)])
            if (2000 <= yi <= 2099) and (1 <= mi <= 12) and (1 <= di <= 31)
            else ""
        )
    except ValueError:
        return ""


def isBenchmarkDir(d):
    with open(os.path.join(d, d + ".status"), "r") as f:
        return str.find(f.readline(), "BENCHMARKS UPDATED") != -1


def findBenchmark(d, t):
    try:
        with open(os.path.join(d, t + ".status"), "r") as f:
            return str.find(f.readline(), "benchmarks updated") != -1
    except:
        return False


def rmDir(d):
    print("  deleting", d)
    shutil.rmtree(d)


if __name__ == "__main__":
    reg_test_gc(sys.argv)
