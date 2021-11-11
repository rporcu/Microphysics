#!/usr/bin/env python3

"""
A simple regression test framework for a AMReX-based code

There are several major sections to this source: the runtime parameter
routines, the test suite routines, and the report generation routines.
They are separated as such in this file.

This test framework understands source based out of the AMReX framework.

"""

from __future__ import print_function

import asyncio
import contextlib
import email
import importlib.util
import logging
import os
import shutil
import smtplib
import subprocess
import sys
import tarfile
import time
import re
import json
from pathlib import Path
from typing import List, Tuple

from rtl2 import params
from rtl2 import test_util
from rtl2 import test_report as report
from rtl2 import test_coverage as coverage
from rtl2.suite import Suite

safe_flags = ["TEST", "USE_CUDA", "USE_ACC", "USE_MPI", "USE_OMP", "DEBUG", "USE_GPU"]


@contextlib.contextmanager
def chdir(path):
    """
    On enter, change directory to specified path.
    On exit, change directory to original.
    """
    this_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(this_dir)


async def main(argv):
    """
    the main test suite driver
    """
    (args, test_list, suite) = setup(argv)
    create_dirs(args, test_list, suite)
    update_time = setup_repos(args, suite)
    if not args.post_only:
        config_time = time.time()
        await suite.cmake_config()
        suite.config_time = time.time() - config_time

    os.chdir(suite.testTopDir)

    # --------------------------------------------------------------------------
    # Get execution times from previous runs
    # --------------------------------------------------------------------------
    runtimes = suite.get_wallclock_history()

    runners = get_runners(args, test_list, suite, runtimes)
    await test_suite(runners, args, suite)
    return finish(args, test_list, suite, update_time, runtimes)


def _active_test_list(test_list):
    return [t.name for t in test_list]


def _check_safety(cs):
    try:
        flag = cs.split("=")[0]
        return flag in safe_flags
    except:
        return False


def check_realclean_safety(compile_strings):
    split_strings = compile_strings.strip().split()
    return all(_check_safety(cs) for cs in split_strings)


def find_build_dirs(tests):
    """given the list of test objects, find the set of UNIQUE build
    directories.  Note if we have the useExtraBuildDir flag set"""

    build_dirs = []
    last_safe = False

    for obj in tests:

        # keep track of the build directory and which source tree it is
        # in (e.g. the extra build dir)

        # first find the list of unique build directories
        dir_pair = (obj.buildDir, obj.extra_build_dir)
        if build_dirs.count(dir_pair) == 0:
            build_dirs.append(dir_pair)

        # re-make all problems that specify an extra compile argument,
        # and the test that comes after, just to make sure that any
        # unique build commands are seen.
        obj.reClean = 1
        if check_realclean_safety(obj.addToCompileString):
            if last_safe:
                obj.reClean = 0
            else:
                last_safe = True

    return build_dirs


def copy_benchmarks(old_full_test_dir, full_web_dir, test_list, bench_dir, log):
    """copy the last plotfile output from each test in test_list
    into the benchmark directory.  Also copy the diffDir, if
    it exists"""
    td = os.getcwd()

    for t in test_list:
        wd = "{}/{}".format(old_full_test_dir, t.name)
        os.chdir(wd)

        if t.compareFile == "" and t.outputFile == "":
            p = t.get_compare_file(output_dir=wd)
        elif not t.outputFile == "":
            if not os.path.exists(t.outputFile):
                p = test_util.get_recent_filename(wd, t.outputFile, ".tgz")
            else:
                p = t.outputFile
        else:
            if not os.path.exists(t.compareFile):
                p = test_util.get_recent_filename(wd, t.compareFile, ".tgz")
            else:
                p = t.compareFile

        if p:
            if p.endswith(".tgz"):
                try:
                    with tarfile.open(name=p, mode="r:gz") as tg:
                        tg.extractall()
                except:
                    log.fail("ERROR extracting tarfile")
                idx = p.rfind(".tgz")
                p = p[:idx]

            store_file = p
            if not t.outputFile == "":
                store_file = "{}_{}".format(t.name, p)

            try:
                shutil.rmtree("{}/{}".format(bench_dir, store_file))
            except:
                pass
            shutil.copytree(p, "{}/{}".format(bench_dir, store_file))

            with open("{}/{}.status".format(full_web_dir, t.name), "w") as cf:
                cf.write("benchmarks updated.  New file:  {}\n".format(store_file))

        else:  # no benchmark exists
            with open("{}/{}.status".format(full_web_dir, t.name), "w") as cf:
                cf.write("benchmarks update failed")

        # is there a diffDir to copy too?
        if not t.diffDir == "":
            diff_dir_bench = "{}/{}_{}".format(bench_dir, t.name, t.diffDir)
            if os.path.isdir(diff_dir_bench):
                shutil.rmtree(diff_dir_bench)
                shutil.copytree(t.diffDir, diff_dir_bench)
            else:
                if os.path.isdir(t.diffDir):
                    try:
                        shutil.copytree(t.diffDir, diff_dir_bench)
                    except IOError:
                        log.warn("file {} not found".format(t.diffDir))
                    else:
                        log.log("new diffDir: {}_{}".format(t.name, t.diffDir))
                else:
                    try:
                        shutil.copy(t.diffDir, diff_dir_bench)
                    except IOError:
                        log.warn("file {} not found".format(t.diffDir))
                    else:
                        log.log("new diffDir: {}_{}".format(t.name, t.diffDir))

        os.chdir(td)


def get_variable_names(suite, plotfile):
    """uses fvarnames to extract the names of variables
    stored in a plotfile"""

    # Run fvarnames
    command = "{} {}".format(suite.tools["fvarnames"], plotfile)
    sout, serr, ierr = test_util.run(command, suite.log)

    if ierr != 0:
        return serr

    # Split on whitespace
    tvars = re.split(r"\s+", sout)[2:-1:2]

    return set(tvars)


def process_comparison_results(stdout, tvars, test):
    """checks the output of fcompare (passed in as stdout)
    to determine whether all relative errors fall within
    the test's tolerance"""

    # Alternative solution - just split on whitespace
    # and iterate through resulting list, attempting
    # to convert the next two items to floats. Assume
    # the current item is a variable if successful.

    # Split on whitespace
    regex = r"\s+"
    words = re.split(regex, stdout)

    indices = filter(lambda i: words[i] in tvars, range(len(words)))

    for i in indices:
        _, _, rel_err = words[i : i + 3]
        if abs(test.tolerance) <= abs(float(rel_err)):
            return False

    return True


def determine_coverage(suite):

    try:
        results = coverage.main(suite.full_test_dir)
    except:
        suite.log.warn("error generating parameter coverage reports, check formatting")
        return

    if all(res is not None for res in results):

        suite.covered_frac = results[0]
        suite.total = results[1]
        suite.covered_nonspecific_frac = results[2]
        suite.total_nonspecific = results[3]

        spec_file = os.path.join(suite.full_test_dir, coverage.SPEC_FILE)
        nonspec_file = os.path.join(suite.full_test_dir, coverage.NONSPEC_FILE)

        shutil.copy(spec_file, suite.full_web_dir)
        shutil.copy(nonspec_file, suite.full_web_dir)


def setup(argv):
    """
    the main test suite driver
    """

    # parse the commandline arguments
    args = test_util.get_args(arg_string=argv)

    # read in the test information
    suite, test_list = params.load_params(args)

    test_list = suite.get_tests_to_run(test_list)

    suite.log.skip()
    suite.log.bold("running tests: ")
    suite.log.indent()
    for obj in test_list:
        suite.log.log(obj.name)
    suite.log.outdent()

    if args.complete_report_from_crash:

        # make sure the web directory from the crash run exists
        suite.full_web_dir = "{}/{}/".format(suite.webTopDir, args.complete_report_from_crash)
        if not os.path.isdir(suite.full_web_dir):
            suite.log.fail("Crash directory does not exist")

        suite.test_dir = args.complete_report_from_crash

        # find all the tests that completed in that web directory
        tests = []
        test_file = ""
        was_benchmark_run = 0
        for sfile in os.listdir(suite.full_web_dir):
            if os.path.isfile(sfile) and sfile.endswith(".status"):
                index = sfile.rfind(".status")
                tests.append(sfile[:index])

                with open(suite.full_web_dir + sfile, "r") as f:
                    for line in f:
                        if line.find("benchmarks updated") > 0:
                            was_benchmark_run = 1

            if os.path.isfile(sfile) and sfile.endswith(".ini"):
                test_file = sfile

        # create the report for this test run
        _num_failed = report.report_this_test_run(
            suite, was_benchmark_run, "recreated report after crash of suite", "", tests, test_file
        )

        # create the suite report
        suite.log.bold("creating suite report...")
        report.report_all_runs(suite, _active_test_list(test_list))
        suite.log.close_log()
        sys.exit("done")

    return (args, test_list, suite)


def create_dirs(args, test_list, suite):

    suite.testTopDir.mkdir(exist_ok=True)
    # --------------------------------------------------------------------------
    # check bench dir and create output directories
    # --------------------------------------------------------------------------
    all_compile = all(t.compileTest == 1 for t in test_list)

    if not all_compile:
        bench_dir = suite.get_bench_dir()

    if args.copy_benchmarks is not None:
        last_run = suite.get_last_run()

    suite.make_test_dirs()

    # Make sure the web dir is valid
    try:
        suite.webTopDir.mkdir(exist_ok=True)
    except OSError:
        suite.log.fail(f"ERROR: unable to create the web directory: {suite.webTopDir}\n")

    if args.copy_benchmarks is not None:
        old_full_test_dir = suite.testTopDir + suite.suiteName + "-tests/" + last_run
        copy_benchmarks(old_full_test_dir, suite.full_web_dir, test_list, bench_dir, suite.log)

        # here, args.copy_benchmarks plays the role of make_benchmarks
        _num_failed = report.report_this_test_run(
            suite,
            args.copy_benchmarks,
            "copy_benchmarks used -- no new tests run",
            "",
            test_list,
            args.input_file[0],
        )
        report.report_all_runs(suite, _active_test_list(test_list))

        sys.exit("done")


def setup_repos(args, suite):
    # --------------------------------------------------------------------------
    # figure out what needs updating and do the git updates, save the
    # current hash / HEAD, and make a ChangeLog
    # --------------------------------------------------------------------------
    now = time.localtime(time.time())
    update_time = time.strftime("%Y-%m-%d %H:%M:%S %Z", now)

    no_update = args.no_update.lower()
    if args.copy_benchmarks is not None:
        no_update = "all"

    # the default is to update everything, unless we specified a hash
    # when constructing the Repo object
    if no_update == "none":
        pass

    elif no_update == "all":
        for k in suite.repos:
            suite.repos[k].update = False

    else:
        nouplist = [k.strip() for k in no_update.split(",")]

        for repo in suite.repos:
            if repo.lower() in nouplist:
                suite.repos[repo].update = False

    os.chdir(suite.testTopDir)

    for k in suite.repos:
        suite.log.skip()
        suite.log.bold("repo: {}".format(suite.repos[k].name))
        suite.log.indent()

        suite.log.outdent()

    # keep track if we are running on any branch that is not the suite
    # default
    branches = [suite.repos[r].get_branch_name() for r in suite.repos]
    if not all(suite.default_branch == b for b in branches):
        suite.log.warn("some git repos are not on the default branch")
        with open("{}/branch.status".format(suite.full_web_dir), "w") as bf:
            bf.write("branch different than suite default")
    return update_time


def get_runners(args, test_list, suite, runtimes):
    return [TestRunner(test, args, test_list, suite, runtimes) for test in test_list]


async def test_suite(runners, args, suite):
    """main loop over tests"""
    for runner in runners:
        runner.check_post_test()

    if not suite.post_only:
        await TestRunner.compile_tests(suite, runners)

        pre_runners = [runner.pre_test() for runner in runners]
        await asyncio.gather(*pre_runners)

        if suite.job_manager == "local":
            await run_local(runners)
        elif suite.job_manager == "slurm":
            await run_slurm(runners, args, suite)
    else:
        for runner in runners:
            runner.test.compile_successful = True

    post_runners = [runner.post_test() for runner in runners]
    await asyncio.gather(*post_runners)


async def run_local(runners):
    the_runners = [runner.run_test() for runner in runners]
    await asyncio.gather(*the_runners)


async def run_slurm(runners, args, suite):
    sb_fname = suite.testDir / "rtl2.conf"
    ntasks, sb_s = srun_script(suite, runners, args)
    assert suite.testDir.is_dir()
    suite.full_test_dir.mkdir(exist_ok=True, parents=True)
    with open(sb_fname, "w") as sb:
        sb.write(sb_s)
    sb_fname.chmod(0o755)
    cmdline = [
        "srun",
        "-W0",
        "-c",
        str(suite.cpus_per_task),
        "-p",
        str(suite.partition),
        "-n",
        str(ntasks),
        "--multi-prog",
        sb_fname.as_posix(),
    ]
    wall_time = time.time()
    subprocess.run(cmdline, check=False)
    wall_time = time.time() - wall_time
    for runner in runners:
        runner.test.wall_time = wall_time


def srun_script(suite: Suite, runners: List["TestRunner"], args: List[str]) -> Tuple[int, str]:
    lines = ["##  srun multiple program configuration file  ##"]
    task_id = 0
    for runner in runners:
        next_id = task_id + runner.test.numprocs
        cmd = runner.test.command(suite, runner.test.base_command(suite, args))
        outdir = runner.test.output_dir
        lines.append(
            f"{task_id}-{next_id}  bash -c '[[ $0 -eq 0 ]] && ( cd {outdir} && {cmd} ) || true' %o"
        )
        task_id = next_id + 1

    ntasks = task_id
    return ntasks, "\n".join(lines)


def toppath():
    return Path(__file__).parent.resolve()


class TestRunner:
    def __init__(self, test, args, test_list, suite, runtimes):
        self.test = test
        self._args = args
        self.test_list = test_list
        self._suite = suite
        self._runtimes = runtimes
        self.test.output_dir = suite.testDir / test.name

    async def pre_test(self):
        suite = self._suite
        test = self.test
        suite.log.outdent()  # just to make sure we have no indentation
        suite.log.skip()
        suite.log.bold(f"working on test: {test.name}")
        suite.log.indent()

        # test.output_dir.mkdir(parents=True, exist_ok=True)

        self.copy_to_rundir()

    async def run_test(self):
        await self.test.run_actual_test(self._suite, self._runtimes, self._args)

    def check_post_test(self):
        prefix = f"{self._suite.Label}." if self._suite.Label else ""
        refdata = self._suite.refdataDir / self.test.name / f"{prefix}runningave.dat"
        if not refdata.parent.is_dir():
            self.test.log.fail(f"Check failed; refdata directory does not exist: {refdata.parent}")
        if not refdata.is_file():
            self.test.log.warn(f"refdata file {refdata} does not exist; creating it...")
            refdata.touch()

    async def post_test(self):
        testname = Path(self.test.output_dir).name
        if self._suite.post_only and not self.test.output_dir.is_dir():
            self.test.log.fail(
                f"Test output dir {self.test.output_dir} does not exist."
                "  You may need to run without --post_only"
            )

        with chdir(self.test.output_dir):
            try:
                testdir = self._suite.sourceDir / self.test.buildDir
                spec = importlib.util.spec_from_file_location(
                    f"run.{testname}.post", testdir / "post.py"
                )
                mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)
                refdata = self._suite.refdata(self.test)
                mod.plot(refdata)
            except Exception:  # pylint: disable=broad-except
                logging.exception("Error post processing %s", testname)

        self.compare()
        self.move_output_to_webdir()
        self.archive_output_and_write_report()

    @staticmethod
    async def compile_tests(suite, runners):
        suite.log.log("building...")

        tests = [runner.test for runner in runners]
        testnames = "_".join(test.name for test in tests)
        output_dir = runners[0].test.output_dir.parent / "cmake_build_tests"
        output_dir.mkdir(parents=True, exist_ok=True)
        coutfile = output_dir / f"{testnames}.build.out"

        prebuild_time = time.time()
        comp_string, rc = await suite.build_tests_cmake(tests=tests, outfile=coutfile)
        build_time = time.time() - prebuild_time

        # copy the build.out into the web directory
        for test in tests:
            (suite.full_web_dir / f"{test.name}.build.out").unlink(missing_ok=True)
            (suite.full_web_dir / f"{test.name}.build.out").symlink_to(coutfile)

        for test in tests:
            test.comp_string = comp_string
            test.build_time = build_time
            test.compile_successful = not rc

            if not test.compile_successful:
                error_msg = "ERROR: compilation failed"
                report.report_single_test(suite, test, tests, failure_msg=error_msg)
                return

            if test.compileTest:
                suite.log.log("creating problem test report ...")
                report.report_single_test(suite, test, tests)

    def copy_to_rundir(self):
        suite = self._suite
        test = self.test
        test_list = self.test_list
        output_dir = test.output_dir
        # ----------------------------------------------------------------------
        # copy the necessary files over to the run directory
        # ----------------------------------------------------------------------
        suite.log.log("copying files to run directory...")

        self.executable = test_util.get_recent_filename(test.bdir(suite), "", ".ex")

        needed_files = []
        if self.executable is not None:
            needed_files.append(self.executable)

        if test.run_as_script:
            needed_files.append(test.run_as_script)

        if test.inputFile:
            needed_files.append(test.inputFile)

        if test.probinFile != "":
            needed_files.append(test.probinFile)

        for auxf in test.auxFiles:
            needed_files.append(auxf)

        # if any copy/move fail, we move onto the next test
        skip_to_next_test = 0
        for nfile in needed_files:
            try:
                dest = output_dir / Path(nfile).name
                if dest.exists():
                    dest.unlink()
                dest.parent.mkdir(exist_ok=True, parents=True)
                os.link(nfile, dest)
            except IOError:
                error_msg = f"ERROR: unable to link file {nfile}"
                report.report_single_test(suite, test, test_list, failure_msg=error_msg)
                skip_to_next_test = 1
                break

        if skip_to_next_test:
            return

        skip_to_next_test = 0
        for lfile in test.linkFiles:
            if not os.path.exists(lfile):
                error_msg = "ERROR: link file {} does not exist".format(lfile)
                report.report_single_test(suite, test, test_list, failure_msg=error_msg)
                skip_to_next_test = 1
                break

            link_source = os.path.abspath(lfile)
            link_name = os.path.join(output_dir, os.path.basename(lfile))
            try:
                os.symlink(link_source, link_name)
            except IOError:
                error_msg = "ERROR: unable to symlink link file: {}".format(lfile)
                report.report_single_test(suite, test, test_list, failure_msg=error_msg)
                skip_to_next_test = 1
                break

        if skip_to_next_test:
            return

    def compare(self):
        suite = self._suite
        test = self.test
        runtimes = self._runtimes
        # ----------------------------------------------------------------------
        # do the comparison
        # ----------------------------------------------------------------------

        suite.log.log(f'Searching "{test.outfile}" for SuccessString: "{test.stSuccessString}"...')

        try:
            with open(test.output_dir / test.outfile, "r") as of:
                # successful comparison is indicated by presence
                # of success string
                for line in of.readlines():
                    if line.find(test.stSuccessString) >= 0:
                        try:
                            wt = line.split()[-1]
                            test.wall_time = float(wt)
                        except ValueError:
                            suite.log.warn(f"Unable to parse walltime from: {line.strip()}")
                        test.compare_successful = True
                        break
        except IOError:
            suite.log.warn("no output file found")

        with open(test.comparison_outfile, "w") as cf:
            if test.compare_successful:
                cf.write("SELF TEST SUCCESSFUL\n")
            else:
                cf.write("SELF TEST FAILED\n")

        # ----------------------------------------------------------------------
        # if the test ran and passed, add its runtime to the dictionary
        # ----------------------------------------------------------------------

        if test.record_runtime(suite):
            test_dict = runtimes.setdefault(test.name, suite.timing_default)
            test_dict["runtimes"].insert(0, test.wall_time)
            test_dict["dates"].insert(0, suite.full_test_dir.name)

    def move_output_to_webdir(self):
        suite = self._suite
        test = self.test
        # ----------------------------------------------------------------------
        # move the output files into the web directory
        # ----------------------------------------------------------------------

        outfile = test.output_dir / test.outfile
        if outfile.is_file():
            (suite.full_web_dir / test.outfile).unlink(missing_ok=True)
            (suite.full_web_dir / test.outfile).symlink_to(outfile)
        else:
            self._suite.log.warn(f"Does not exist:  {outfile}")
        if os.path.isfile(test.errfile):
            (suite.full_web_dir / test.errfile.name).unlink(missing_ok=True)
            (suite.full_web_dir / test.errfile.name).symlink_to(test.errfile)
            test.has_stderr = True
        if test.doComparison:
            shutil.copy(test.comparison_outfile, suite.full_web_dir)
        analysis = Path(f"{test.name}.analysis.out")
        if analysis.is_file():
            shutil.copy(analysis, suite.full_web_dir)

        inputs = test.inputFile
        if inputs:
            shutil.copy(inputs, suite.full_web_dir / f"{test.name}.{inputs.name}")

        if test.probinFile != "":
            shutil.copy(test.probinFile, suite.full_web_dir / f"{test.name}.{test.probinFile}")

        for af in test.auxFiles:
            shutil.copy(af, suite.full_web_dir / f"{test.name}.{af.name}")

        if test.png_file is not None:
            try:
                shutil.copy(test.png_file, suite.full_web_dir)
            except IOError:
                # visualization was not successful.  Reset image
                test.png_file = None

        if not test.analysisRoutine == "":
            try:
                shutil.copy(test.analysisOutputImage, suite.full_web_dir)
            except IOError:
                suite.log.warn("unable to copy analysis image")
                # analysis was not successful.  Reset the output image
                test.analysisOutputImage = ""

        # were any Backtrace files output (indicating a crash)
        suite.copy_backtrace(test)

        out_pngfile = test.output_dir / "output.png"
        web_pngfile = (suite.full_web_dir / test.name).with_suffix(".png")
        if out_pngfile.is_file():
            shutil.copy(out_pngfile, web_pngfile)

    def archive_output_and_write_report(self):
        suite = self._suite
        test = self.test
        test_list = self.test_list
        args = self._args
        output_dir = test.output_dir
        output_file = ""
        # ----------------------------------------------------------------------
        # archive (or delete) the output
        # ----------------------------------------------------------------------
        suite.log.log("archiving the output...")
        for pfile in output_dir.iterdir():

            if pfile.is_dir() and re.match(f"{test.name}.*_(plt|chk)[0-9]+", pfile.name):

                if suite.purge_output == 1 and not pfile == output_file:

                    # delete the plt/chk file
                    try:
                        shutil.rmtree(pfile)
                    except:
                        suite.log.warn(f"unable to remove {pfile}")

                else:
                    # tar it up
                    try:
                        with tarfile.open("{}.tgz".format(pfile), "w:gz") as tar:
                            tar.add(str(pfile))

                    except:
                        suite.log.warn(f"unable to tar output file {pfile}")

                    else:
                        try:
                            shutil.rmtree(pfile)
                        except OSError:
                            suite.log.warn(f"unable to remove {pfile}")

        # ----------------------------------------------------------------------
        # write the report for this test
        # ----------------------------------------------------------------------
        if args.make_benchmarks is None:
            suite.log.log("creating problem test report ...")
            report.report_single_test(suite, test, test_list)


def finish(args, test_list, suite, update_time, runtimes):

    # --------------------------------------------------------------------------
    # jsonify and save runtimes
    # --------------------------------------------------------------------------
    file_path = suite.get_wallclock_file()
    with open(file_path, "w") as json_file:
        json.dump(runtimes, json_file, indent=4)

    # --------------------------------------------------------------------------
    # parameter coverage
    # --------------------------------------------------------------------------
    if suite.reportCoverage:
        determine_coverage(suite)

    # --------------------------------------------------------------------------
    # write the report for this instance of the test suite
    # --------------------------------------------------------------------------
    suite.log.outdent()
    suite.log.skip()
    suite.log.bold("creating new test report...")
    num_failed = report.report_this_test_run(
        suite, args.make_benchmarks, args.note, update_time, test_list, args.input_file[0]
    )

    # make sure that all of the files in the web directory are world readable
    for file_ in os.listdir(suite.full_web_dir):
        current_file = suite.full_web_dir / file_
        if current_file.is_file():
            os.chmod(current_file, 0o644)

    # For temporary run, return now without creating suite report.
    if args.do_temp_run:
        suite.delete_tempdirs()
        return num_failed

    # store an output file in the web directory that can be parsed easily by
    # external program
    name = "source"
    branch = "HEAD"
    with open("{}/suite.{}.status".format(suite.webTopDir, branch.replace("/", "_")), "w") as f:
        f.write(
            "{}; num failed: {}; source hash: {}".format(
                suite.repos[name].name, num_failed, suite.repos[name].hash_current
            )
        )

    # --------------------------------------------------------------------------
    # generate the master report for all test instances
    # --------------------------------------------------------------------------
    suite.log.skip()
    suite.log.bold("creating suite report...")
    report.report_all_runs(suite, _active_test_list(test_list))

    # delete any temporary directories
    suite.delete_tempdirs()

    def email_developers():
        msg = email.message_from_string(suite.emailBody)
        msg["From"] = suite.emailFrom
        msg["To"] = ",".join(suite.emailTo)
        msg["Subject"] = suite.emailSubject

        server = smtplib.SMTP("localhost")
        server.sendmail(suite.emailFrom, suite.emailTo, msg.as_string())
        server.quit()

    if num_failed > 0 and suite.sendEmailWhenFail and not args.send_no_email:
        suite.log.skip()
        suite.log.bold("sending email...")
        email_developers()

    return num_failed


if __name__ == "__main__":
    todo = main(sys.argv[1:])

    loop = asyncio.get_event_loop()
    n = loop.run_until_complete(todo)

    sys.exit(n)
