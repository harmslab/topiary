#!/usr/bin/env python
"""
Crawl through source and test directories to identify tests that are missing.
"""

import os, re, sys, string

def code_loader(filename,is_test):
    """
    Get list of all functions in filename.

    Parameters
    ----------
    filename : str
        file to load
    is_test : bool
        whether or not this is test (look for def test_) or no (look for def )

    Returns
    -------
    def_lines : list
        list of tuples. first element is function, second is whether function
        does anything (or just has a "pass")
    """

    if is_test:
        function_finder = re.compile("def test_.*?\\(.*?")
    else:
        function_finder = re.compile("def .*?\\(.*?")

    class_check = re.compile("class .*?:")
    in_class_check = re.compile("def .*?\\(self")

    def_lines = []
    current_function = None
    in_docstring = False
    with open(filename) as f:
        for line in f:

            # Is this a class?
            if class_check.search(line):
                def_lines.append((line.strip(),True))
                continue

            # Is this a function?
            if function_finder.search(line):

                # Skip functions inside other functions by looking for large
                # indent. Use class check to keep methods to classes
                if line.startswith("    "):
                    if not in_class_check.search(line):
                        continue

                current_function = line.strip()
                continue

            # Now strip line
            line = line.strip()

            if current_function:

                # Blank line or comment
                if line.split("#")[0].strip() == "":
                    continue

                # Checking for whether we're in a docstring
                if line.startswith("\"\"\""):

                    if in_docstring:
                        in_docstring = False
                    else:
                        in_docstring = True
                    continue

                # If we're currently in a docstring, keep going.
                if in_docstring:
                    continue

                # Okay, now have a non-docstring, non-blank, non-comment line
                if line.strip() == "pass":
                    def_lines.append((current_function,False))
                    current_function = None
                else:
                    def_lines.append((current_function,True))
                    current_function = None

    return def_lines


def crawler(some_directory,is_test):
    """
    Crawl through code in some directory and return all functions seen.
    """

    cache_pattern = re.compile("__pycache__")

    functions = {}
    for dirpath, dirs, files in os.walk(some_directory):

        if cache_pattern.search(dirpath):
            continue

        for f in files:
            if f[-3:] != ".py":
                continue

            filename = os.path.join(dirpath,f)
            function_lines = code_loader(filename,is_test)
            functions[filename] = function_lines

    return functions

def clean_dir_name(dir_name):
    if dir_name.startswith(os.path.sep):
        dir_name = dir_name.strip(os.path.sep)
        dir_name = f"{os.path.sep}{dir_name}"
    else:
        dir_name = dir_name.strip(os.path.sep)
    return dir_name

def completeness_crawler(code_dir,test_dir,bin_dir=None):
    """
    Compare expected and found tests.
    """

    code_dir = clean_dir_name(code_dir)
    test_dir = clean_dir_name(test_dir)
    if bin_dir is not None:
        bin_dir = clean_dir_name(bin_dir)

    # These are going to have an extra test inserted for uniqueness. For example,
    # A.util.a and B.util.b will have tests A.test_A_util.test_a and
    # B.test_B_util.test_b rather than A.test_util.test_a and B.test_util.test_b
    repeated_modules = ["util","model","tree"]

    # python code in tests that is not really tests
    skip_test_modules = ["conftest","completeness_crawler"]

    # Get functions from both codebase and tests
    code_functions = crawler(code_dir,is_test=False)
    test_functions = crawler(test_dir,is_test=False)

    all_functions = []
    dummy_functions = []

    # For each file...
    for k in code_functions:

        base = re.sub(os.path.sep,".",k)[:-3]

        # For each function...
        class_name = None
        for elements in code_functions[k]:

            # Get function name and whether it has more than just "pass"
            f, is_real = elements

            # If this is a "class" line, indicate we are in a class. We're
            # going to want to test base.class_name...
            if f.startswith("class"):
                class_name = f.split(" ")[1].split("(")[0].strip(":")
                all_functions.append(f"{base}.{class_name}")
                continue

            # Get the function name
            fcn = f.split(" ")[1].split("(")[0]

            # If this is an __init__ function, skip it.
            if fcn == "__init__":
                continue

            # If we have a class defined...
            if class_name is not None:

                # Record as method
                if re.search("\\(self",f):
                    fcn = f"{class_name}.{fcn}"

                # Otherwise, not a method, we've moved out of the class
                else:
                    class_name = None

            # resolved funciton name
            fcn = f"{base}.{fcn}"

            # If this is a real function, append it to all functions. Ignore
            # fake functions; they don't need to be tested.
            if is_real:
                all_functions.append(fcn)
            else:
                dummy_functions.append(fcn)

    # Now figure out expected test names for each function
    expected_tests = []
    for a in all_functions:

        # topiary.blah.yo.function -> tests.blah.test_yo.test_function
        chunks = a.split(".")
        chunks[0] = "tests"

        # Is previous chunk a class? If so, combine into single test with name
        # test_ClassName_method. Trim chunks back so repeated_modules bit works
        # below.
        if re.search("[" + string.ascii_uppercase + "]",chunks[-2]):
            chunks[-2] = f"test_{chunks[-2]}_{chunks[-1]}"
            chunks.pop(-1)
        else:
            chunks[-1] = f"test_{chunks[-1]}"

        # for repeated modules, add submodule into test nomenclature
        if chunks[-2] in repeated_modules:
            chunks[-2] = f"test_{chunks[-3]}_{chunks[-2]}"
        else:
            chunks[-2] = f"test_{chunks[-2]}"

        # Record expected tests
        expected_tests.append(".".join(chunks))

    extra_functions = []
    real_tests = []
    fake_tests = []

    # Now go through the test modules
    for k in test_functions:

        base = re.sub(os.path.sep,".",k)[:-3]

        # Skip things like conftest
        if base.split(".")[-1] in skip_test_modules:
            continue

        # Go through the test functions
        for elements in test_functions[k]:

            f, is_real = elements

            fcn = f.strip().split(" ")[1].split("(")[0]

            # If it starts with "test_" it's a test function
            if fcn.startswith("test_"):

                # Record whether or not it's a real test
                if is_real:
                    real_tests.append(f"{base}.{fcn}")
                else:
                    fake_tests.append(f"{base}.{fcn}")

            # Otherwise, it's just some other function to worry about
            else:
                extra_functions.append(f"{base}.{fcn}")

    # Get any expected script tests (say bin/run-it)
    if bin_dir is not None:

        scripts = os.listdir(bin_dir)

        for s in scripts:
            test = f"{test_dir}.{bin_dir}.test_{s}.test_main"
            expected_tests.append(test)

    real_tests = set(real_tests)
    expected_tests = set(expected_tests)
    fake_tests = set(fake_tests)

    good_tests = expected_tests.intersection(real_tests)
    all_tests_seen = real_tests.union(fake_tests)

    missing_tests = expected_tests - all_tests_seen
    extra_tests = list(all_tests_seen - expected_tests)

    integration_tests = []
    truly_extra_tests = []
    for e in extra_tests:
        if re.search("test_integrated",e):
            integration_tests.append(e)
        else:
            truly_extra_tests.append(e)

    extra_fake_tests = fake_tests - expected_tests
    missing_fake_tests = fake_tests.intersection(missing_tests)

    good_tests = list(good_tests)
    good_tests.sort()
    for g in good_tests:
        print("FOUND",g)

    integration_tests = list(integration_tests)
    integration_tests.sort()
    for i in integration_tests:
        print("FOUND & INTEGRATION",i)

    fake_tests = list(fake_tests)
    fake_tests.sort()
    for f in fake_tests:
        print("FAKE",f)

    missing_tests = list(missing_tests)
    missing_tests.sort()
    for m in missing_tests:
        print("MISSING",m)

    truly_extra_tests = list(truly_extra_tests)
    truly_extra_tests.sort()
    for e in truly_extra_tests:
        print("FOUND & EXTRA",e)

    return missing_tests

def write_test_file(filename,functions):

    chunks = filename[:-3].split(os.path.sep)[1:]
    chunks[-1] = "_".join(chunks[-1].split("_")[1:])
    to_import = ".".join(chunks)

    import_functions = []
    for test_function in functions:
        import_functions.append(re.sub("test_","",test_function))

    if not os.path.exists(filename):

        out = []
        out.append("import pytest")
        out.append("import topiary")
        out.append("")
        for f in import_functions:
            out.append(f"from topiary.{to_import} import {f}")
        out.append("")

        f = open(filename,'w')
        f.write("\n".join(out))
        f.close()

    f = open(filename,'a')
    for fcn in functions:
        f.write(f"def {fcn}():\n\n    pass\n\n")
    f.close()


def template_tests(missing_tests):

    missing_tests = list(missing_tests)
    missing_tests.sort()

    to_write = {}
    for m in missing_tests:
        chunks = m.split(".")

        if re.search("[" + string.ascii_uppercase + "]",chunks[-2]):
            chunks[-2] = f"test_{chunks[-2]}_{chunks[-1]}"
            chunks.pop(-1)

        function = chunks[-1]
        filename = f"{os.path.sep.join(chunks[:-1])}.py"

        try:
            to_write[filename].append(function)
        except KeyError:
            to_write[filename] = [function]

    # for filename in to_write:
    #     write_test_file(filename,to_write[filename])


def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    try:
        code_dir = argv[0]
        test_dir = argv[1]
    except IndexError:
        err = "Incorrect arguments. Should be `test_crawler.py code_dir test_dir [bin_dir]`\n"
        raise ValueError(err)

    try:
        bin_dir = argv[2]
    except IndexError:
        bin_dir = None

    missing_tests = completeness_crawler(code_dir,test_dir,bin_dir)

    template_tests(missing_tests)


if __name__ == "__main__":
    main()
