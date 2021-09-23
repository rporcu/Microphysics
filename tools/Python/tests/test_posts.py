from pathlib import Path
import contextlib
import importlib
import os

import pytest

test_files = {
    "uio_p_g_1.dat": """#  Time   p_g  vol
0 1.691 6.21629e-10
0.00114004 5.35665 6.21629e-10
0.00237973 7.20442 6.21629e-10
0.003586 8.54199 6.21629e-10""",
    "uio_p_g_2.dat": """#  Time   p_g  vol
0 1.63535 6.21629e-10
0.00114004 5.16973 6.21629e-10
0.00237973 6.94743 6.21629e-10
0.003586 8.22839 6.21629e-10""",
    "uio_vel_p_0.dat": """#  Time   NP  U  V  W  KE
0 7986 0 0 0 0
0.00114004 7986 -0.00670169 -2.6395e-06 -3.66679e-06 1.18674e-14
0.00237973 7986 -0.0123184 -9.30656e-06 -8.55365e-06 4.05456e-14
0.003586 7986 -0.016514 -1.07292e-05 -1.46242e-05 7.38519e-14""",
}


@pytest.mark.xfail(reason="FIXME")
@pytest.mark.parametrize(
    "testname",
    [
        ("fluidbed_cyl"),
        ("fluidbed_periodic"),
        ("fluidbed_sq"),
        ("hopper"),
        ("riser_cyl"),
        ("riser_sq"),
    ],
)
def test_post(testname, tmp_path):
    for test_file, content in test_files.items():
        with open(tmp_path / test_file, "w") as test_f:
            test_f.write(content)

    with chdir(tmp_path):

        mod = importlib.import_module(f"run.{testname}.post")
        new_png = mod.plot()
        assert new_png == Path("output.png")


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
