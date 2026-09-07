import pytest
import os
import tempfile
import socket
import subprocess

from netguard import NetworkAccessBlocked
from netguard import is_local_address

# Markers that mean "the caller explicitly opted into external access". A test
# behind any of these is only collected when its flag is passed, so its use of
# the network is a stated cost, not a hidden one -- the guard leaves it alone.
#
# run_generax and run_raxml are in this list because several of them really do
# reach the Open Tree of Life API: setting up a GeneRax run annotates the
# species tree from OTL. That is surprising, but it is not *hidden*, and
# blocking it would just break --run-generax in CI.
_OPT_IN_MARKERS = frozenset(["network",
                             "run_ncbi_server",
                             "run_blast",
                             "run_generax",
                             "run_raxml"])

# pytester lets tests/test_harness.py run throwaway pytest sessions to verify
# the autouse fixtures below actually contain cross-test contamination.
pytest_plugins = ["pytester"]

def pytest_addoption(parser):
    """
    Add options to the pytest command line parser.
    """

    parser.addoption("--run-generax",
                     action="store_true",
                     default=False,
                     help="Run tests involving generax")

    parser.addoption("--run-raxml",
                     action="store_true",
                     default=False,
                     help="Run tests involving raxml")

    parser.addoption("--run-blast",
                     action="store_true",
                     default=False,
                     help="Run tests involving blast")

    parser.addoption("--run-ncbi-server",
                     action="store_true",
                     default=False,
                     help="Run tests that access the NCBI server")

    parser.addoption("--run-network",
                     action="store_true",
                     default=False,
                     help="Run tests that require live internet access")

    parser.addoption("--fd-report",
                     action="store_true",
                     default=False,
                     help="Track open file descriptors per test and print the "
                          "biggest leakers at the end of the session")

def pytest_collection_modifyitems(config, items):
    """
    Look for run_generax, run_raxml, run_blast, and run_ncbi_server decorators. 
    Modify test collection based on 1) pytest command line arguments and 
    2) operating system.
    """

    # Look for --run-generax argument. Skip test if this is not specified.
    if not config.getoption("--run-generax"):
        skipper = pytest.mark.skip(reason="Only run when --run-generax is given")
        for item in items:
            if "run_generax" in item.keywords:
                item.add_marker(skipper)

    # Look for --run-raxml argument. Skip test if this is not specified.
    if not config.getoption("--run-raxml"):
        skipper = pytest.mark.skip(reason="Only run when --run-raxml is given")
        for item in items:
            if "run_raxml" in item.keywords:
                item.add_marker(skipper)

    # Look for --run-blast argument. Skip test if this is not specified.
    if not config.getoption("--run-blast"):
        skipper = pytest.mark.skip(reason="Only run when --run-blast is given")
        for item in items:
            if "run_blast" in item.keywords:
                item.add_marker(skipper)

    # Look for --run-ncbi-server argument. Skip test if this is not specified.
    if not config.getoption("--run-ncbi-server"):
        skipper = pytest.mark.skip(reason="Only run when --run-ncbi-server is given")
        for item in items:
            if "run_ncbi_server" in item.keywords:
                item.add_marker(skipper)

    # Look for --run-network argument. Skip test if this is not specified.
    if not config.getoption("--run-network"):
        skipper = pytest.mark.skip(reason="Only run when --run-network is given")
        for item in items:
            if "network" in item.keywords:
                item.add_marker(skipper)

    # Assign every test exactly one tier marker: unit, smoke, or integration.
    #
    # `integration` is derived, never hand-written: a test that needs an
    # external binary or a live service is already behind an opt-in flag, and
    # deriving the tier from that flag means the two can never disagree.
    #
    # `smoke` is written in the source, because "uses real data / real file
    # I/O and takes a moment" is a judgement, not something collection can see.
    #
    # `unit` is the default. A new test is hermetic and fast until someone says
    # otherwise -- and the block_subprocess fixture below enforces that, so an
    # unlabelled test that shells out fails instead of quietly slowing the
    # fast suite down.
    for item in items:

        if len(set(item.keywords) & _OPT_IN_MARKERS) > 0:
            item.add_marker(pytest.mark.integration)
        elif "smoke" in item.keywords:
            pass
        else:
            item.add_marker(pytest.mark.unit)

    # If this is a windows box, skip any test with run_generax or run_raxml
    # decorators.
    if os.name == "nt":
        disallowed_dec = set(["run_generax","run_raxml"])
        skipper = pytest.mark.skip(reason="cannot run on windows")
        for item in items:
            if len(set(item.keywords).intersection(disallowed_dec)) > 0:
                item.add_marker(skipper)

@pytest.fixture(autouse=True)
def restore_cwd():
    """
    Put the process back in its starting directory after every test, whether
    the test passed, failed, or raised.

    Many tests do `current_dir = os.getcwd(); os.chdir(tmpdir); ...` and then
    restore on the last line of the body. That restore never runs if anything
    above it fails, so a single failure strands the interpreter's cwd inside a
    temp directory that pytest is about to delete -- and every later test in
    the session inherits it. That is why the suite produced path-dependent
    failures that did not reproduce when tests were run individually.

    This fixture makes the restore unconditional, so no test can leak its
    working directory into the next one.
    """

    start_dir = os.getcwd()

    try:
        yield
    finally:
        try:
            os.chdir(start_dir)
        except OSError:
            # The starting directory itself went away (a test deleted it).
            # Anything valid beats staying in a directory that no longer
            # exists, since os.getcwd() would then raise for everyone after us.
            os.chdir(tempfile.gettempdir())


def _open_fd_count():
    """
    Number of file descriptors this process currently holds, or None on a
    platform that does not expose them this way.
    """

    for fd_dir in ("/proc/self/fd", "/dev/fd"):
        try:
            return len(os.listdir(fd_dir))
        except OSError:
            continue

    return None


# test nodeid -> file descriptors opened and never closed by that test
_FD_GROWTH = {}

# Running tally, written as the session goes so the data survives a crash.
_FD_LOG = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "reports", "fd-running-total.txt")
_FD_LOG = os.path.abspath(_FD_LOG)


@pytest.fixture(autouse=True)
def track_file_descriptors(request):
    """
    With --fd-report, record how many file descriptors each test leaves behind.

    Running the whole suite can exhaust the process file-descriptor limit while
    running the same tests individually does not. That points at a test (or a
    piece of topiary) that opens descriptors and never closes them, so the cost
    only shows up once hundreds of tests share one process.

    This does not fail anything -- it just measures, so the next full run
    (`pytest tests --run-raxml --run-generax --run-blast --fd-report`) names
    the culprit instead of leaving it to guesswork.
    """

    if not request.config.getoption("--fd-report"):
        yield
        return

    before = _open_fd_count()
    try:
        yield
    finally:
        after = _open_fd_count()
        if before is not None and after is not None:
            if after > before:
                _FD_GROWTH[request.node.nodeid] = after - before

            # Write as we go. Descriptor exhaustion kills the session before
            # pytest_terminal_summary can run, so a report that only prints at
            # the end is exactly the report you cannot get when you need it.
            try:
                with open(_FD_LOG, "a") as f:
                    f.write(f"{after}\t{after - before:+d}\t{request.node.nodeid}\n")
            except OSError:
                pass


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """
    Say what pytest is exiting with, and print the file-descriptor report if
    one was requested.

    The exit status is reported whenever it is not success. A green summary
    followed by a nonzero exit status is indistinguishable, in a CI log, from a
    wrapper (coverage.py, the shell) inventing a failure of its own -- and a
    nightly that reported "452 passed" and then failed the build with exit 1
    cost an afternoon to tell those two apart. One line settles it: if the line
    is absent, pytest returned 0 and the failure came from outside pytest.
    """

    if exitstatus != 0:
        terminalreporter.write_line(f"pytest exit status: {int(exitstatus)}")

    if not config.getoption("--fd-report"):
        return

    if _open_fd_count() is None:
        terminalreporter.write_line(
            "file descriptors could not be counted on this platform")
        return

    terminalreporter.section("file descriptor growth")

    if len(_FD_GROWTH) == 0:
        terminalreporter.write_line("no test leaked a file descriptor")
        return

    total = sum(_FD_GROWTH.values())
    terminalreporter.write_line(
        f"{len(_FD_GROWTH)} tests leaked {total} descriptors in total")

    worst = sorted(_FD_GROWTH.items(), key=lambda kv: kv[1], reverse=True)
    for nodeid, grew in worst[:25]:
        terminalreporter.write_line(f"  +{grew:<4} {nodeid}")


@pytest.fixture(autouse=True)
def block_network(request):
    """
    Stop tests from silently reaching the live internet.

    19 tests in the default run used to call out to the Open Tree of Life API,
    ftp.ncbi.nlm.nih.gov, and a CDN -- with no marker saying so. That made a
    green run mean "the network was up" rather than "the code is correct", and
    it is the most likely source of the intermittent CI failures.

    The guard applies to tests that are *not* behind an opt-in flag -- that is,
    the default suite, which is where all 19 hidden dependencies lived. Tests
    that genuinely need the internet without needing an external binary are
    marked `network` and are skipped unless --run-network is given. Tests
    already gated behind --run-generax, --run-raxml, --run-blast or
    --run-ncbi-server are left alone (see _OPT_IN_MARKERS): the caller asked
    for them, so their network use is declared rather than hidden.

    The effect is that a newly-introduced hidden dependency in a supposedly
    hermetic test fails loudly and by name, instead of becoming the next
    intermittent CI failure.
    """

    if len(set(request.keywords) & _OPT_IN_MARKERS) > 0:
        yield
        return

    real_connect = socket.socket.connect
    real_connect_ex = socket.socket.connect_ex

    def guarded_connect(self, address, *args, **kwargs):
        if is_local_address(self, address):
            return real_connect(self, address, *args, **kwargs)
        raise NetworkAccessBlocked(
            f"{request.node.nodeid} tried to connect to {address!r}.\n"
            f"If this test really needs the internet, mark it "
            f"@pytest.mark.network (it will then be skipped unless "
            f"--run-network is given). Otherwise, mock the call.")

    def guarded_connect_ex(self, address, *args, **kwargs):
        if is_local_address(self, address):
            return real_connect_ex(self, address, *args, **kwargs)
        raise NetworkAccessBlocked(
            f"{request.node.nodeid} tried to connect to {address!r}.")

    socket.socket.connect = guarded_connect
    socket.socket.connect_ex = guarded_connect_ex

    try:
        yield
    finally:
        socket.socket.connect = real_connect
        socket.socket.connect_ex = real_connect_ex


class SubprocessBlocked(Exception):
    """
    A test in the `unit` tier tried to launch an external process.
    """

    pass


@pytest.fixture(autouse=True)
def block_subprocess(request):
    """
    Keep the `unit` tier honest by making it impossible to shell out.

    A tier is only useful if it means something. `unit` claims "hermetic and
    fast"; without enforcement that claim decays the first time someone adds a
    test that quietly runs muscle. Blocking subprocess creation turns that from
    a slow fast-suite into an immediate, named failure.

    A test that legitimately needs to run something external belongs in `smoke`
    (mark it `@pytest.mark.smoke`) or, if it needs a real external binary,
    behind one of the --run-* flags.
    """

    if "unit" not in request.keywords:
        yield
        return

    real_popen = subprocess.Popen

    def guarded_popen(*args, **kwargs):
        cmd = args[0] if args else kwargs.get("args")
        raise SubprocessBlocked(
            f"{request.node.nodeid} is in the `unit` tier but tried to run "
            f"{cmd!r}.\nMark it @pytest.mark.smoke, or put it behind a "
            f"--run-* flag if it needs an external binary.")

    subprocess.Popen = guarded_popen
    try:
        yield
    finally:
        subprocess.Popen = real_popen
