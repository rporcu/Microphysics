import os
from pathlib import Path
from typing import Optional, Tuple

import rtl2.test_coverage as coverage
from rtl2.suite import Suite
from rtl2.test_util import git_commit


CSS_CONTENTS = r"""
body {font-family: "Arial", san-serif;}

h1 {font-family: "Tahoma","Arial", sans-serif;
    color: #333333;}

h3 {display: inline;}

h3.passed {text-decoration: none; display: inline;
           color: black; background-color: lime; padding: 2px;}

a.passed:link {color: black; text-decoration: none;}
a.passed:visited {color: black; text-decoration: none;}
a.passed:hover {color: #ee00ee; text-decoration: underline;}

a.passed-slowly:link {color: black; text-decoration: none;}
a.passed-slowly:visited {color: black; text-decoration: none;}
a.passed-slowly:hover {color: #ee00ee; text-decoration: underline;}

h3.failed {text-decoration: none; display: inline;
           color: yellow; background-color: red; padding: 2px;}

a.failed:link {color: yellow; text-decoration: none;}
a.failed:visited {color: yellow; text-decoration: none;}
a.failed:hover {color: #00ffff; text-decoration: underline;}

a.compfailed:link {color: yellow; text-decoration: none;}
a.compfailed:visited {color: yellow; text-decoration: none;}
a.compfailed:hover {color: #00ffff; text-decoration: underline;}

a.crashed:link {color: yellow; text-decoration: none;}
a.crashed:visited {color: yellow; text-decoration: none;}
a.crashed:hover {color: #00ffff; text-decoration: underline;}

h3.benchmade {text-decoration: none; display: inline;
              color: black; background-color: orange; padding: 2px;}

a.benchmade:link {color: black; text-decoration: none;}
a.benchmade:visited {color: black; text-decoration: none;}
a.benchmade:hover {color: #00ffff; text-decoration: underline;}

span.nobreak {white-space: nowrap;}
span.mild-success {color: green;}
span.mild-failure {color: red;}

a.main:link {color: yellow; text-decoration: none;}
a.main:visited {color: yellow; text-decoration: none;}
a.main:hover {color: #00ffff; text-decoration: underline;}

td {border-width: 0px;
    padding: 5px;
    background-color: white;
    vertical-align: middle;}

td.passed {background-color: lime; opacity: 0.8;}
td.passed-slowly {background-color: yellow; opacity: 0.8;}
td.failed {background-color: red; color: yellow; opacity: 0.8;}
td.compfailed {background-color: purple; color: yellow; opacity: 0.8;}
td.crashed {background-color: black; color: yellow; opacity: 0.8;}
td.benchmade {background-color: orange; opacity: 0.8;}
td.date {background-color: #666666; color: white; opacity: 0.8; font-weight: bold;}

.maintable tr:hover {background-color: blue;}


table {border-collapse: separate;
       border-spacing: 2px;
       margin-left: auto;
       margin-right: auto;
       border-width: 1px;
       border-color: gray;
       border-style: solid;
       box-shadow: 10px 10px 5px #888888;}

table.head {border-collapse: separate;
       border-spacing: 0px;
       margin-left: auto;
       margin-right: auto;
       border-width: 0px;
       border-style: solid;
       box-shadow: none;}

/* http://blog.petermares.com/2010/10/27/vertical-text-in-html-table-headers-for-webkitmozilla-browsers-without-using-images/ */

div.verticaltext {text-align: center;
                  vertical-align: middle;
                  width: 20px;
                  margin: 0px;
                  padding: 0px;
                  padding-left: 3px;
                  padding-right: 3px;
                  padding-top: 10px;
                  white-space: nowrap;
                  -webkit-transform: rotate(-90deg);
                  -moz-transform: rotate(-90deg);}

#summary th {background-color: grey;
    color: yellow;
    text-align: center;
    height: 2em;
    padding-bottom: 3px;
    padding-left: 5px;
    padding-right: 5px;}


#summary td {background: transparent;}

#summary tr:nth-child(even) {background: #dddddd;}
#summary tr:nth-child(odd) {background: #eeeeee;}

#summary tr.special {background: #ccccff;}
#summary td.highlight {color: red;}

#summary td.passed {background-color: lime; }
#summary td.passed-slowly {background-color: yellow; }
#summary td.failed {background-color: red; color: yellow;}
#summary td.benchmade {background-color: orange;}
#summary td.compfailed {background-color: purple; color: yellow;}
#summary td.crashed {background-color: black; color: yellow;}

div.small {font-size: 75%;}

th {background-color: grey;
    color: yellow;
    text-align: center;
    vertical-align: bottom;
    height: @TABLEHEIGHT@;
    padding-bottom: 3px;
    padding-left: 5px;
    padding-right: 5px;}

li {padding-top: 0.5em;}

ul li {color: blue;
       font-weight: bold;}
ul li ul li {color: black;
             font-weight: normal;}

ul li h3 {border: 1px solid black;}

#compare td {font-family: "Lucida Console", Monaco, monospace;
             font-size: 80%;}

#box {  width: 900px;
  margin: 0 auto;
  padding: 1em;
  background: #ffffff;
}

.alignright {
   text-align: right;
}

"""

HTML_HEADER = r"""
<HTML>
<HEAD>
<TITLE>@TESTDIR@ / @TESTNAME@</TITLE>
<META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=ISO-8859-1">
<LINK REL="stylesheet" TYPE="text/css" HREF="tests.css">
</HEAD>
<BODY>
<div id="box">
"""

MAIN_HEADER = r"""
<HTML>
<HEAD>
<TITLE>@TITLE@</TITLE>
<META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=ISO-8859-1">
<LINK REL="stylesheet" TYPE="text/css" HREF="tests.css">
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.5.0/css/font-awesome.min.css">
</HEAD>
<BODY>
<!--GOUPLINK-->
<CENTER><H1>@TITLE@</H1></CENTER>
<CENTER><H2>@SUBTITLE@</H2></CENTER>
<P><TABLE class='maintable'>
<CENTER>
  <td align=center class="benchmade"><h3>Benchmark Updated</h3></td>
  <td align=center class="failed"><h3>Comparison Failed</h3></td>
  <td align=center class="compfailed"><h3>Compilation Failed</h3></td>
  <td align=center class="crashed"><h3>Crashed</h3></td>
  <td align=center class="passed"><h3>Passed</h3></td>
  <td align=center class="passed-slowly"><h3>Performance Drop</h3></td>
</CENTER>
</TABLE>
"""


def create_css(table_height=16):
    """write the css file for the webpages"""

    css = CSS_CONTENTS.replace("@TABLEHEIGHT@", "{}em".format(table_height))

    with open("tests.css", "w") as cf:
        cf.write(css)


class HTMLList:
    """a simple class for managing nested HTML lists"""

    def __init__(self, of=None):
        # items will hold tuples: (indent, string), where indent
        # specifies how deeply nested we are
        self.list_items = []
        self.current_indent = 0
        self.of = of

    def item(self, content):
        # add an item to the list
        self.list_items.append((self.current_indent, content))

    def indent(self):
        # indent (nest a new list)
        self.current_indent += 1

    def outdent(self):
        # close the current nest level
        self.current_indent -= 1

    def write_list(self):
        # output the list to the outfile, of, specified at creation
        self.of.write("<ul>\n")
        current_indent = -1
        for i, c in self.list_items:
            if current_indent == -1:
                current_indent = i
            else:
                if i < current_indent:
                    self.of.write("</li></ul></li>\n")
                elif i > current_indent:
                    self.of.write("<ul>\n")
                else:
                    self.of.write("</li>\n")

            current_indent = i
            self.of.write("<li>{}\n".format(c))

        # finish current item
        self.of.write("</li>")

        # finish nesting
        for _ in range(0, current_indent):
            self.of.write("</ul></li>\n")

        self.of.write("</ul>\n")


class HTMLTable:
    """a simple class for creating an HTML table"""

    def __init__(self, out_file, columns=1, divs=None):
        """create the table object.  Here divs is the name of
        any HTML div(s) we want to wrap the table with"""

        self.hf = out_file
        self.columns = columns
        if divs is not None:
            self.divs = list(divs)
        else:
            self.divs = None

    def start_table(self):
        if self.divs is not None:
            for d in self.divs:
                self.hf.write("<div id='{}'>\n".format(d))
        self.hf.write("<p><table>\n")

    def header(self, header_list):
        """write the table header"""
        n = len(header_list)
        line = "<tr>" + n * "<th>{}</th>" + "</tr>\n"
        self.hf.write(line.format(*header_list))

    def print_single_row(self, row):
        self.hf.write(f"<tr class='special'><td colspan={self.columns}>{row}</td></tr>\n")

    def print_row(self, row_list, highlight=False):
        """row_list are the individual table elements.  Note that if
        a list item is a tuple, then the first element is assumed to
        be the cell data and the second element is an html tag that
        goes in the <td >, e.g. to set the class or colspan"""

        n = len(row_list)
        if highlight:
            line = "<tr>" + n * "<td class='highlight'>{}</td>" + "</tr>\n"
        else:
            line = "<tr>"
            for d in row_list:
                if isinstance(d, tuple):
                    line += "<td {}>{}</td>".format(d[1], d[0])
                else:
                    line += "<td>{}</td>".format(d)
            line += "</tr>\n"
        self.hf.write(line.format(*row_list))

    def end_table(self):
        self.hf.write("</table>\n")
        if self.divs is not None:
            for _ in range(len(self.divs)):
                self.hf.write("</div>\n")


def get_particle_compare_command(diff_lines):
    for line in diff_lines:
        if line.find("particle_compare") > 0:
            return line
    return ""


def report_single_test(suite, test, tests, failure_msg=None):
    """generate a single problem's test result page.  If
    failure_msg is set to a string, then it is assumed
    that the test did not complete.  The string will
    be reported on the test page as the error."""

    # for navigation
    tnames = [t.name for t in tests]
    current_index = tnames.index(test.name)

    if failure_msg is not None:
        suite.log.testfail("aborting test")
        suite.log.testfail(failure_msg)

    current_dir = os.getcwd()
    os.chdir(suite.full_web_dir)

    # we stored compilation success in the test object
    compile_successful = test.compile_successful

    analysis_successful = True
    if test.analysisRoutine != "":
        analysis_successful = test.analysis_successful

    # we store comparison success in the test object but also read
    # in the comparison report for displaying
    if failure_msg is None:
        if not test.compileTest:
            compare_successful = test.compare_successful

            if test.doComparison:
                compare_file = test.comparison_outfile
                try:
                    with open(compare_file, "r") as cf:
                        diff_lines = cf.readlines()
                except IOError:
                    suite.log.warn("WARNING: no comparison file found")
                    diff_lines = [""]

            # last check: did we produce any backtrace files?
            if test.crashed:
                compare_successful = False

        # write out the status file for this problem, with either
        # PASSED, PASSED SLOWLY, COMPILE FAILED, or FAILED
        status_file = "{}.status".format(test.name)
        with open(status_file, "w") as sf:
            if not compile_successful:
                sf.write("COMPILE FAILED\n")
                suite.log.testfail(f"{test.name} COMPILE FAILED")
            elif test.compileTest or (compare_successful and analysis_successful):
                string = "PASSED\n"
                if test.check_performance:
                    meets_threshold, _, _ = test.measure_performance()
                    if not (meets_threshold is None or meets_threshold):
                        string = "PASSED SLOWLY\n"
                sf.write(string)
                suite.log.success(f"{test.name} PASSED")
            elif test.crashed:
                sf.write("CRASHED\n")
                suite.log.testfail(f"{test.name} CRASHED (backtraces produced)")
            else:
                sf.write("FAILED\n")
                suite.log.testfail(f"{test.name} FAILED.  For details see {test.outfile}")

    else:
        # we came in already admitting we failed...
        if not test.compile_successful:
            msg = "COMPILE FAILED"
        else:
            msg = "FAILED"

        status_file = "{}.status".format(test.name)
        with open(status_file, "w") as sf:
            sf.write("{}\n".format(msg))
        suite.log.testfail("{} {}".format(test.name, msg))

    # --------------------------------------------------------------------------
    # generate the HTML page for this test
    # --------------------------------------------------------------------------

    # write the css file
    create_css()

    html_file = "{}.html".format(test.name)
    hf = open(html_file, "w")

    new_head = HTML_HEADER

    # arrows for previous and next test
    new_head += r"""<table style="width: 100%" class="head"><br><tr>"""
    if current_index > 0:
        new_head += r"""<td><< <a href="{}.html">previous test</td>""".format(
            tests[current_index - 1].name
        )
    else:
        new_head += r"""<td>&nbsp;</td>"""

    if current_index < len(tests) - 1:
        new_head += r"""<td class="alignright"><a href="{}.html">next test >></td>""".format(
            tests[current_index + 1].name
        )
    else:
        new_head += r"""<td>&nbsp;</td>"""

    new_head += r"</tr></table>" + "\n"

    new_head += r"""<center><h1><a href="index.html">@TESTDIR@</a> / @TESTNAME@</h1></center>"""

    new_head = new_head.replace("@TESTDIR@", suite.full_test_dir.as_posix())
    new_head = new_head.replace("@TESTNAME@", test.name)

    hf.write(new_head)

    ll = HTMLList(of=hf)

    if failure_msg is not None:
        ll.item("Test error: ")
        ll.indent()

        ll.item('<h3 class="failed">Failed</h3>')
        ll.item("{}".format(failure_msg))

        ll.outdent()

    # build summary
    ll.item("Build/Test information:")
    ll.indent()

    ll.item("Build directory: {}".format(test.buildDir))

    if not test.extra_build_dir == "":
        ll.indent()
        ll.item("in {}".format(suite.repos[test.extra_build_dir].dir))
        ll.outdent()

    if not test.compileTest:

        if test.debug:
            ll.item("Debug test")

        if test.acc:
            ll.item("OpenACC test")

        if test.useMPI or test.useOMP:
            ll.item("Parallel run")
            ll.indent()
            if test.useMPI:
                ll.item("MPI numprocs = {}".format(test.numprocs))
            if test.useOMP:
                ll.item("OpenMP numthreads = {}".format(test.numthreads))
            ll.outdent()

        if test.restartTest:

            ll.item("Restart test")
            ll.indent()
            ll.item(
                f"Job was run as normal and then restarted from checkpoint # {test.restartFileNum},"
                " and the two final outputs were compared"
            )
            ll.outdent()

        ll.item("Files:")
        ll.indent()

        if test.inputFile:
            ll.item(f"input file: <a href='{test.inputFile}'>{test.name}.{test.inputFile.name}</a>")

        if test.probinFile != "":
            ll.item(
                'probin file: <a href="{}.{}">{}</a>'.format(
                    test.name, test.probinFile, test.probinFile
                )
            )

        for i, afile in enumerate(test.auxFiles):
            # sometimes the auxFile was in a subdirectory under the
            # build directory.
            root_file = os.path.basename(afile)
            ll.item(
                'auxillary file {}: <a href="{}.{}">{}</a>'.format(
                    i + 1, test.name, root_file, afile
                )
            )

        ll.outdent()

        ll.item("Dimensionality: {}".format(test.dim))

    ll.outdent()  # end of build information

    # compilation summary
    ll.item("Compilation:")
    ll.indent()

    if compile_successful:
        ll.item('<h3 class="passed">Successful</h3>')
    else:
        ll.item('<h3 class="failed">Failed</h3>')

    ll.item("Configuration time: {:.3f} s".format(test.config_time))
    ll.item("Configuration command:<br><tt>{}</tt>".format(test.config_cmd))
    ll.item('<a href="{}.cmake.log">configuration output</a>'.format(suite.suiteName))
    ll.item("Compilation time: {:.3f} s".format(test.build_time))
    ll.item("Compilation command:<br><tt>{}</tt>".format(test.comp_string))
    ll.item('<a href="{}.build.out">build output</a>'.format(test.name))

    ll.outdent()

    if not test.compileTest:

        # execution summary
        ll.item("Execution:")
        ll.indent()
        ll.item("Execution time: {:.3f} s".format(test.wall_time))

        if test.check_performance:

            meets_threshold, percentage, compare_str = test.measure_performance()

            if meets_threshold is not None:

                if meets_threshold:
                    style = "mild-success"
                else:
                    style = "mild-failure"

                ll.item("{} run average: {:.3f} s".format(test.runs_to_average, test.past_average))
                ll.item(
                    'Relative performance: <span class="{}">{:.1f}% {}</span>'.format(
                        style, percentage, compare_str
                    )
                )

        ll.item("Execution command:<br><tt>{}</tt>".format(test.run_command))
        ll.item('<a href="{}.run.out">execution output</a>'.format(test.name))
        if test.has_stderr:
            ll.item('<a href="{}.err.out">execution stderr</a>'.format(test.name))
        if test.has_jobinfo:
            ll.item('<a href="{}.job_info">job_info</a>'.format(test.name))
        ll.outdent()

        # were there backtrace files?
        if test.crashed:
            ll.item("Backtraces:")
            ll.indent()
            for bt in test.backtrace:
                ll.item(f'<a href="{bt}">{bt}</a>')
            ll.outdent()

        # comparison summary
        if failure_msg is None:
            ll.item("Comparison: ")
            ll.indent()

            if compare_successful:
                ll.item('<h3 class="passed">Successful</h3>')
            else:
                ll.item('<h3 class="failed">Failed</h3>')
            ll.outdent()

        if test.analysisRoutine != "":
            ll.item("Analysis: ")
            ll.indent()

            if test.analysis_successful:
                ll.item('<h3 class="passed">Successful</h3>')
            else:
                ll.item('<h3 class="failed">Failed</h3>')

            ll.item('<a href="{}.analysis.out">execution output</a>'.format(test.name))
            ll.outdent()

        ll.item('<img src="{}.png" alt="Comparison vs historical results"/>'.format(test.name))

    ll.write_list()

    if (not test.compileTest) and test.doComparison and failure_msg is None:

        # parse the compare output and make an HTML table
        ht = HTMLTable(hf, columns=3, divs=["summary", "compare"])
        in_diff_region = False

        box_error = False
        grid_error = False
        variables_error = False
        no_bench_error = False
        particle_counts_differ_error = False

        pcomp_line = get_particle_compare_command(diff_lines)

        for line in diff_lines:
            if "number of boxes do not match" in line:
                box_error = True
                break

            if "grids do not match" in line:
                grid_error = True
                break

            if "number of variables do not match" in line:
                variables_error = True

            if "no corresponding benchmark found" in line:
                no_bench_error = True
                break

            if "Particle data headers do not agree" in line:
                particle_counts_differ_error = True
                break

            if in_diff_region:
                # diff region
                hf.write(line)
                continue

            if line.find("fcompare") > 1:
                hf.write("<tt>" + line + "</tt>\n")
                if pcomp_line:
                    hf.write("<tt>" + pcomp_line + "</tt>\n")

                ht.start_table()
                continue

            if line.strip().startswith("diff "):
                # this catches the start of a plain text diff --
                # we need the space here to not match variables
                # that start with diff
                ht.end_table()
                hf.write("<pre>\n")

                hf.write(line)
                in_diff_region = True
                continue

            if line.strip().startswith("level "):
                ht.print_single_row(line.strip())
                continue

            if line.strip().startswith("-----"):
                continue

            if line.strip().startswith("<<<"):
                ht.print_single_row(line.strip().replace("<", "&lt;").replace(">", "&gt;"))
                continue

            fields = [q.strip() for q in line.split("  ") if not q == ""]

            if fields:
                if fields[0].startswith("variable"):
                    ht.header(fields)
                    continue

                if len(fields) == 2:
                    if "NaN present" in line:
                        ht.print_row([fields[0], (fields[1], "colspan='2'")])
                        continue
                    if "variable not present" in line:
                        ht.print_row([fields[0], (fields[1], "colspan='2'")])
                        continue
                    ht.header([" "] + fields)
                    continue

                if len(fields) == 1:
                    continue

                abs_err = float(fields[1])
                rel_err = float(fields[2])
                if abs(rel_err) > 1.0e-6:
                    ht.print_row([fields[0], abs_err, rel_err], highlight=True)
                else:
                    ht.print_row([fields[0], abs_err, rel_err])

        if in_diff_region:
            hf.write("</pre>\n")
        else:
            ht.end_table()

        if box_error:
            hf.write("<p>number of boxes do not match</p>\n")

        if grid_error:
            hf.write("<p>grids do not match</p>\n")

        if no_bench_error:
            hf.write("<p>no corresponding benchmark found</p>\n")

        if variables_error:
            hf.write("<p>variables differ in files</p>\n")

        if particle_counts_differ_error:
            hf.write("<p>number of particles differ in files</p>\n")

    if (not test.compileTest) and failure_msg is None:
        # show any visualizations
        if test.doVis:
            if test.png_file is not None:
                hf.write("<P>&nbsp;\n")
                hf.write("<P><IMG SRC='{}' BORDER=0>".format(test.png_file))

        # show any analysis
        if not test.analysisOutputImage == "":
            hf.write("<P>&nbsp;\n")
            hf.write("<P><IMG SRC='%s' BORDER=0>" % (test.analysisOutputImage))

    # close
    hf.write("</div></body>\n")
    hf.write("</html>\n")

    hf.close()

    # switch back to the original directory
    os.chdir(current_dir)


def report_this_test_run(suite, make_benchmarks, note, _update_time, test_list, test_file):
    """generate the master page for a single run of the test suite"""

    # get the current directory
    current_dir = os.getcwd()

    # switch to the web directory and open the report file
    os.chdir(suite.full_web_dir)

    try:
        build_time = max([q.build_time for q in test_list])
    except:
        build_time = -1

    try:
        wall_time = max([q.wall_time for q in test_list])
    except:
        wall_time = -1

    # keep track of the number of tests that passed and the number that failed
    num_failed = 0
    num_passed = 0

    # --------------------------------------------------------------------------
    # generate the HTML page for this run of the test suite
    # --------------------------------------------------------------------------

    # always create the css (in case it changes)
    create_css()

    with open("LABEL", "w") as label_f:
        label_f.write(suite.Label)

    # create the master web page
    hf = open("index.html", "w")

    new_head = (
        HTML_HEADER + r"""<CENTER><H1><A HREF="../">@TESTDIR@</A> / @TESTNAME@</H1></CENTER>"""
    )

    new_head = new_head.replace("@TESTDIR@", suite.suiteName)
    new_head = new_head.replace("@TESTNAME@", suite.full_test_dir.as_posix())

    hf.write(new_head)

    if not note == "":
        hf.write('<p><b>Test run note:</b><br><font color="gray">%s</font>\n' % (note))

    if make_benchmarks is not None:
        hf.write(
            '<p><b>Benchmarks updated</b><br>comment: <font color="gray">{}</font>\n'.format(
                make_benchmarks
            )
        )

    hf.write('<p><b>test input parameter file:</b> <A HREF="%s">%s</A>\n' % (test_file, test_file))

    if build_time > 0:
        hf.write("<p><b>build time for all tests:</b> {} s\n".format(build_time))

    if wall_time > 0:
        hf.write("<p><b>wall clock time for all tests:</b> {} s\n".format(wall_time))

    mfix_branch, mfix_hash = git_commit(suite.sourceDir)
    amrex_branch, amrex_hash = git_commit(suite.sourceDir / "subprojects" / "amrex")
    csgeb_branch, csgeb_hash = git_commit(suite.sourceDir / "subprojects" / "csg-eb")
    hydro_branch, hydro_hash = git_commit(suite.sourceDir / "subprojects" / "AMReX-Hydro")
    amrex_url = f"https://github.com/AMReX-Codes/amrex/commit/{amrex_hash}"
    hydro_url = f"https://github.com/AMReX-Codes/AMReX-Hydro/commit/{hydro_hash}"
    csgeb_url = f"https://mfix.netl.doe.gov/gitlab/exa/csg-eb/-/commits/{csgeb_hash}"
    mfix_url = f"https://mfix.netl.doe.gov/gitlab/exa/mfix/-/commits/{mfix_hash}"
    hf.write("<ul>")
    hf.write(
        f"<li><b>AMReX</b>; <b>branch:</b> {amrex_branch}; "
        f"<b>hash:</b> <a href='{amrex_url}'>{amrex_hash}</a></li>"
    )
    hf.write(
        f"<li><b>AMReX-Hydro</b>; <b>branch:</b> {hydro_branch}; "
        f"<b>hash:</b> <a href='{hydro_url}'>{hydro_hash}</a></li>"
    )
    hf.write(
        f"<li><b>CSG-EB</b>; <b>branch:</b> {csgeb_branch}; "
        f"<b>hash:</b> <a href='{csgeb_url}'>{csgeb_hash}</a></li>"
    )
    hf.write(
        f"<li><b>MFiX</b>; <b>branch:</b> {mfix_branch}; "
        f"<b>hash:</b> <a href='{mfix_url}'>{mfix_hash}</a></li>"
    )
    hf.write("</ul>")

    hf.write("<p>&nbsp;\n")

    # summary table
    if make_benchmarks is None:
        special_cols = []
        cols = (
            [
                "test name",
                "dim",
                "compare plotfile",
                "# levels",
                "MPI procs",
                "OMP threads",
                "OpenACC",
                "debug",
                "compile",
                "restart",
            ]
            + special_cols
            + ["build time", "wall time", "result"]
        )
        ht = HTMLTable(hf, columns=len(cols), divs=["summary"])
        ht.start_table()
        ht.header(cols)

    else:
        ht = HTMLTable(hf, columns=3, divs=["summary"])
        ht.start_table()
        ht.header(["test name", "result", "comment"])

    # loop over the tests and add a line for each
    for test in test_list:
        if make_benchmarks is None:
            status_file = Path("%s.status" % (test.name)).resolve()
            if not status_file.is_file():
                suite.log.fail(f"Unable to find {status_file.as_posix()}")
            status = None
            with open(status_file, "r") as sf:
                for line in sf:
                    if line.find("PASSED") >= 0:
                        status = "passed"
                        td_class = "passed-slowly" if "SLOWLY" in line else "passed"
                        num_passed += 1
                    elif line.find("COMPILE FAILED") >= 0:
                        status = "compile fail"
                        td_class = "compfailed"
                        num_failed += 1
                    elif line.find("CRASHED") >= 0:
                        status = "crashed"
                        td_class = "crashed"
                        num_failed += 1
                    elif line.find("FAILED") >= 0:
                        status = "failed"
                        td_class = "failed"
                        num_failed += 1

                    if status is not None:
                        break

            row_info = []
            row_info.append('<a href="{}.html">{}</a>'.format(test.name, test.name))
            row_info.append(test.dim)
            row_info.append("<div class='small'>{}</div>".format(test.compare_file_used))

            if test.nlevels is not None:
                row_info.append(test.nlevels)
            else:
                row_info.append("")

            if test.useMPI:
                row_info.append("&check; ({})".format(test.numprocs))
            else:
                row_info.append("")

            # OMP ?
            if test.useOMP:
                row_info.append("&check; ({})".format(test.numthreads))
            else:
                row_info.append("")

            # OpenACC ?
            if test.acc:
                row_info.append("&check;")
            else:
                row_info.append("")

            # debug ?
            if test.debug:
                row_info.append("&check;")
            else:
                row_info.append("")

            # compile ?
            if test.compileTest:
                row_info.append("&check;")
            else:
                row_info.append("")

            # restart ?
            if test.restartTest:
                row_info.append("&check;")
            else:
                row_info.append("")

            # build time
            row_info.append("{:.3f}&nbsp;s".format(test.build_time))

            # wallclock time
            row_info.append("{:.3f}&nbsp;s".format(test.wall_time))

            # result
            row_info.append((status.upper(), "class='{}'".format(td_class)))

            ht.print_row(row_info)

        else:
            if test.restartTest:
                continue
            if test.compileTest:
                continue
            if test.selfTest:
                continue

            # the benchmark was updated -- find the name of the new benchmark file
            benchStatusFile = "%s.status" % (test.name)

            bench_file = "none"

            with open(benchStatusFile, "r") as bf:
                for line in bf:
                    index = line.find("file:")
                    if index >= 0:
                        bench_file = line[index + 5 :]
                        break

            row_info = []
            row_info.append("{}".format(test.name))
            if bench_file != "none":
                row_info.append(("BENCHMARK UPDATED", "class='benchmade'"))
                row_info.append("new benchmark file is {}".format(bench_file))
            else:
                row_info.append(("BENCHMARK NOT UPDATED", "class='failed'"))
                row_info.append("compilation or execution failed")

            ht.print_row(row_info)

    ht.end_table()

    # Test coverage
    if suite.reportCoverage:
        report_coverage(hf, suite)

    # close
    hf.write("</div></body>\n")
    hf.write("</html>\n")
    hf.close()

    # --------------------------------------------------------------------------
    # write out a status file for all the tests
    # --------------------------------------------------------------------------

    status_file = suite.full_test_dir.with_suffix(".status")
    with open(status_file, "w") as sf:

        if make_benchmarks is None:
            if num_failed == 0:
                sf.write("ALL PASSED\n")
            elif num_failed > 0 and num_passed > 0:
                sf.write("SOME FAILED\n")
            else:
                sf.write("ALL FAILED\n")

        else:
            sf.write("BENCHMARKS UPDATED\n")

    # switch back to the original directory
    os.chdir(current_dir)

    return num_failed


def report_coverage(html_file, suite):

    tvars = (
        suite.covered_frac,
        suite.total,
        suite.covered_nonspecific_frac,
        suite.total_nonspecific,
    )
    if not all(tvars):
        return

    cols = ["coverage type", "coverage %", "# covered", "# uncovered"]
    ht = HTMLTable(html_file, len(cols), divs=["summary"])

    ht.start_table()
    ht.header(cols)

    # Overall coverage
    row_info = []
    row_info.append('<a href="{}">{}</a>'.format(coverage.SPEC_FILE, "overall"))
    row_info.append("{:.2f}%".format(100 * suite.covered_frac))
    covered = int(round(suite.total * suite.covered_frac))
    uncovered = suite.total - covered
    row_info.append("{}".format(covered))
    row_info.append("{}".format(uncovered))
    ht.print_row(row_info)

    # Nonspecific-only coverage
    row_info = []
    row_info.append('<a href="{}">{}</a>'.format(coverage.NONSPEC_FILE, "nonspecific only"))
    row_info.append("{:.2f}%".format(100 * suite.covered_nonspecific_frac))
    covered = int(round(suite.total_nonspecific * suite.covered_nonspecific_frac))
    uncovered = suite.total_nonspecific - covered
    row_info.append("{}".format(covered))
    row_info.append("{}".format(uncovered))
    ht.print_row(row_info)

    ht.end_table()


def report_all_runs(suite, active_test_list, max_per_page=50):

    table_height = min(max(suite.lenTestName, 4), 18)

    os.chdir(suite.webTopDir)

    create_css(table_height=table_height)

    valid_dirs, all_tests = suite.get_run_history(active_test_list)

    if suite.do_timings_plots:
        suite.make_timing_plots(all_tests=all_tests)

    # how many pages are we going to spread this over?
    npages = int(len(valid_dirs) / max_per_page) + 1

    for n in range(npages):

        # --------------------------------------------------------------------------
        # generate the HTML
        # --------------------------------------------------------------------------
        title = "%s regression tests" % (suite.suiteName)

        if n == 0:
            hf = open("index.html", "w")
        else:
            hf = open("index{}.html".format(n), "w")

        lvalid_dirs = valid_dirs[
            n * max_per_page : min((n + 1) * max_per_page, len(valid_dirs) - 1)
        ]
        lvalid_dirs = valid_dirs

        header = MAIN_HEADER.replace("@TITLE@", title).replace("@SUBTITLE@", suite.sub_title)

        if suite.goUpLink:
            header2 = header.replace("<!--GOUPLINK-->", '<a href="../">GO UP</a>')
            hf.write(header2)
        else:
            hf.write(header)

        hf.write("<P><TABLE class='maintable'>\n")

        # write out the header
        hf.write("<TR><TH ALIGN=CENTER>date</TH>\n")
        for test in all_tests:
            hf.write("<TH><div class='verticaltext'>%s</div></TH>\n" % (test))

        hf.write("</TR>\n")

        if suite.do_timings_plots:
            hf.write("<tr><td class='date'>plots</td>")
            for t in all_tests:
                plot_file = "{}-timings.{}".format(t, suite.plot_ext)
                if os.path.isfile(plot_file):
                    hf.write(
                        f'<TD ALIGN=CENTER title="{t} timings plot"><H3>'
                        f'<a href="{plot_file}"><i class="fa fa-line-chart"></i></a></H3></TD>\n'
                    )
                else:
                    hf.write("<TD ALIGN=CENTER><H3>&nbsp;</H3></TD>\n")

            hf.write("</TR>\n")

        # loop over all the test runs
        for lvalid_dir in lvalid_dirs:
            tdir = lvalid_dir.relative_to(suite.webTopDir)
            label = ""
            label_path = lvalid_dir / "LABEL"
            if label_path.is_file():
                with open(label_path) as label_f:
                    label = label_f.read().strip()
            # write out the directory (date)
            branch_mark = ""
            hf.write(
                "<TR><TD class='date'><SPAN CLASS='nobreak'>"
                f"<A class='main' HREF=\"{tdir}/index.html\">{tdir}&nbsp;{label}&nbsp;</A>"
                f"{branch_mark}</SPAN></TD>\n"
            )

            for test in all_tests:
                # look to see if the current test was part of this suite run
                status_file = (lvalid_dir / test).with_suffix(".status")
                status, emoji = get_result(suite, status_file)

                # write out this test's status
                if status is None:
                    hf.write("<td>&nbsp;</td>\n")
                elif status == "benchmade":
                    hf.write(
                        '<td align=center title="{}" class="{}"><h3>U</h3></td>\n'.format(
                            test, status
                        )
                    )
                else:
                    hf.write(
                        f'<td align=center title="{test}" class="{status}">'
                        f'<h3><a href="{tdir}/{test}.html" class="{status}">{emoji}</a></h3></td>\n'
                    )

            hf.write("</TR>\n\n")

        hf.write("</TABLE>\n")

        if n != npages - 1:
            hf.write('<p><a href="index{}.html">older tests</a>'.format(n + 1))

        # close
        hf.write("</BODY>\n")
        hf.write("</HTML>\n")

        hf.close()


def get_result(suite: Suite, status_file: Path) -> Tuple[Optional[str], str]:
    if not status_file.is_file():
        suite.log.warn(f"Missing status file {status_file}")
        return (None, "!&nbsp;")

    with open(status_file, "r") as sf:
        for line in sf:
            result = get_line_result(line)
            if result[0] is not None:
                return result
    return (None, "!&nbsp;")


def get_line_result(line: str) -> Tuple[Optional[str], str]:
    return (
        (("passed-slowly", ":]") if "SLOWLY" in line else ("passed", ":)"))
        if "PASSED" in line
        else ("compfailed", ":(")
        if "COMPILE FAILED" in line
        else ("crashed", "xx")
        if "CRASHED" in line
        else ("failed", "!&nbsp;")
        if "FAILED" in line
        else ("benchmade", "U")
        if "benchmarks updated" in line
        else (None, "!&nbsp;")
    )
