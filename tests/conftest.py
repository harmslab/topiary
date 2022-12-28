import pytest
import pandas as pd
import os
import glob
import inspect
import json
from html.parser import HTMLParser

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

def pytest_collection_modifyitems(config, items):
    """
    Look for run_generax and run_raxml decorators. Modify test collection based
    on 1) pytest command line arguments and 3) operating system.
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

    # If this is a windows box, skip any test with run_generax or run_raxml
    # decorators.
    if os.name == "nt":
        disallowed_dec = set(["run_generax","run_raxml"])
        skipper = pytest.mark.skip(reason="cannot run on windows")
        for item in items:
            if len(set(item.keywords).intersection(disallowed_dec)) > 0:
                item.add_marker(skipper)


def get_files(base_dir):
    """
    Traverse base_dir and return a dictionary that keys all files and some
    rudimentary *.ext expressions to absolute paths to those files. They keys
    will be things like "some_dir/test0/rocket.txt" mapping to
    "c:\some_dir\life\base_dir\some_dir\test0\rocket.txt". The idea is to have
    easy-to-read cross-platform keys within unit tests.

    Classes of keys:

        + some_dir/test0/rocket.txt maps to a file (str)
        + some_dir/test0/ maps to the test0 directory itself (str)
        + some_dir/test0/*.txt maps to all .txt (list)
        + some_dir/test0/* maps to all files or directories in the directory
          (list)

    Note that base_dir is *not* included in the keys. All are relative to that
    directory by :code:`os.path.basename(__file__)/{base_dir}`.

    Parameters
    ----------
    base_dir : str
        base directory for search. should be relative to test file location.

    Returns
    -------
    output : dict
        dictionary keying string paths to absolute paths
    """

    containing_dir = os.path.dirname(os.path.realpath(__file__))
    starting_dir = os.path.abspath(os.path.join(containing_dir,base_dir))

    base_length = len(starting_dir.split(os.sep))

    # Traverse starting_dir
    output = {}
    for root, dirs, files in os.walk(starting_dir):

        # path relative to base_dir as a list
        this_path = root.split(os.sep)[base_length:]

        # Build paths to specific files
        local_files = []
        for file in files:
            local_files.append(os.path.join(root,file))
            new_key = this_path[:]
            new_key.append(file)
            output["/".join(new_key)] = local_files[-1]

        # Build paths to patterns of file types
        patterns = {}
        ext = list(set([f.split(".")[-1] for f in local_files]))
        for e in ext:
            new_key = this_path[:]
            new_key.append(f"*.{e}")
            output["/".join(new_key)] = glob.glob(os.path.join(root,f"*.{e}"))

        # Build path to all files in this directory
        new_key = this_path[:]
        new_key.append("*")
        output["/".join(new_key)] = glob.glob(os.path.join(root,f"*"))

        # Build paths to directories in this directory
        for this_dir in dirs:
            new_key = this_path[:]
            new_key.append(this_dir)
            # dir without terminating /
            output["/".join(new_key)] = os.path.join(root,this_dir)

            # dir with terminating /
            new_key.append("")
            output["/".join(new_key)] = os.path.join(root,this_dir)

    return output

def get_public_param_defaults(public_function,private_function):
    """
    Get the defaults for a public function in the API that then passes
    a subset of its arguments directly into a private function. This allows
    us to test the defaults on the public api without having to hardcode
    them into the tests.

    Parameters
    ----------
    public_function : function
        a publically facing function that has kwargs with default values
    private_function : function
        an internal function that has some subset of those kwargs, with or
        without defaults

    Returns
    -------
    kwargs : dict
        kwargs to pass into the private function containing defaults taken from
        the public function.
    """

    public_param = inspect.signature(public_function).parameters
    private_param = inspect.signature(private_function).parameters

    kwargs = {}
    for p in private_param:
        try:
            default = public_param[p].default
        except KeyError:
            # If here, private function has a private argument
            continue

        if default is not inspect._empty:
            kwargs[p] = public_param[p].default

    return kwargs

class HTMLValidator(HTMLParser):
    """
    Class for validating html.
    """

    DOES_NOT_STACK = ["br","img","meta","link"]

    def __init__(self,*args,**kwargs):
        
        self._stack = []
        self._tag_dict = {}
        self._counter = 0

        super(HTMLValidator, self).__init__(*args, **kwargs)
    
    def handle_starttag(self, tag, attrs):
        """
        Start tag.
        """
    
        # If something like "img" or "meta" that does not have a closing tag,
        # don't stick in the stack. 
        if not tag in self.DOES_NOT_STACK:

            # Record tag in stack
            self._stack.append(tag)
        
        # Create dictionary from tag attributes
        attr_dict = {}
        for a in attrs:
            
            key = a[0]
            values = a[1].strip().split()
            attr_dict[key] = values
        
        # Create new entry in tag_dict if needed
        if tag not in self._tag_dict:
            self._tag_dict[tag] = []
        
        # tag_dict is keyed to tag, then records positions, attributes, and data
        # for each element seen. data starts as None. 
        self._tag_dict[tag].append([self._counter,attr_dict,None])
        self._counter += 1
    
        
    def handle_endtag(self, tag):
        """
        End tag.
        """
        
        # Delete last tag from stack, only if it matches. If the last tag does
        # not match, there's mangled html. 
        if self._stack[-1] == tag:
            self._stack = self._stack[:-1]
        else:
            err = "mangled html\n"
            err += f"observed </{tag}> but expected </{self._stack[-1]}>\n"
            raise RuntimeError(err)

    def handle_data(self, data):
        """
        Data
        """
        
        # Skip empty data
        if data.strip() == "":
            return None

        # If data comes in, append it to the last tag
        last_tag = self._stack[-1]
        self._tag_dict[last_tag][-1][-1] = data
    

    @property
    def is_valid(self):
        """
        Properly nested html will have nothing in the stack.
        """

        return len(self._stack) == 0

    @property
    def tag_dict(self):
        """
        Dictionary of tags. Keys are tags. Values are order seen in feed
        (integer), dictionary of this tag attributes, and (possibly) data.
        """
        
        return self._tag_dict

@pytest.fixture(scope="module")
def ncbi_lines():
    """
    A list of generic ncbi lines.
    """
    dir = os.path.dirname(os.path.realpath(__file__))
    filename = os.path.join(dir,"data","ncbi_lines","ncbi-lines-file.txt")

    lines = []
    with open(filename) as f:
        for line in f:
            lines.append(line.strip())

    filename = os.path.join(dir,"data","ncbi_lines","ncbi-lines-parsed.json")
    f = open(filename)
    parsed_lines =json.load(f)
    f.close()

    return (lines, parsed_lines)


@pytest.fixture(scope="module")
def dataframe_good_files():
    """
    A list of files that should be properly parsable by topiary.
    """

    return get_files(os.path.join("data","supported-df-formats"))


@pytest.fixture(scope="module")
def recip_blast_hit_dfs():
    """
    Load saved hit_dfs output.
    """

    def _load(prefix):

        dir = os.path.dirname(os.path.realpath(__file__))
        files = glob.glob(os.path.join(dir,
                                       "data",
                                       "recip_blast",
                                       f"{prefix}_hit_dfs",
                                       "*.csv"))
        files.sort()
        for i in range(len(files)):
            files[i] = pd.read_csv(files[i])

        return files

    return {"ncbi":_load("ncbi"),"local":_load("local")}


@pytest.fixture(scope="module")
def test_dataframes():
    """
    A dictionary holding dataframes of varying sorts of badness to test parser.
    """

    dir = os.path.dirname(os.path.realpath(__file__))
    search_string = os.path.join(dir,"data","test_dataframes","*.csv")

    df_dict = {}
    for g in glob.glob(search_string):
        key = ".".join(os.path.split(g)[1].split(".")[:-1])
        df_dict[key] = pd.read_csv(g)

    return df_dict


@pytest.fixture(scope="module")
def programs():
    """
    Dictionary holding paths pointing to programs to run.
    """

    dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.abspath(os.path.join(dir,"data","programs"))
    files = os.listdir(base_dir)

    out_dict = {}
    for f in files:
        out_dict[f] = os.path.join(base_dir,f)

    return out_dict

@pytest.fixture(scope="module")
def xml():
    """
    Dictionary holding paths pointing to ncbi style xml files.
    """

    dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.abspath(os.path.join(dir,"data","xml"))
    files = os.listdir(base_dir)

    out_dict = {}
    for f in files:
        out_dict[f] = os.path.join(base_dir,f)

    return out_dict

@pytest.fixture(scope="module")
def user_xml_files():
    """
    Dictionary holding paths pointing to ncbi style xml files.
    """

    dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.abspath(os.path.join(dir,"data","xml","user-xml-files"))
    json_file = os.path.join(base_dir,"index.json")

    with open(json_file) as f:
        file_info = json.load(f)

    final_dict = {}
    for f in file_info:
        final_dict[os.path.join(base_dir,f)] = file_info[f]

    return final_dict




@pytest.fixture(scope="module")
def ncbi_blast_server_output():
    """
    These csv files are output from topiary.ncbi.blast.ncbi._thread_manager
    """

    dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.abspath(os.path.join(dir,"data","ncbi-blast-server-output"))

    files = glob.glob(os.path.join(base_dir,"*.csv"))
    files.sort()

    all_hits = []
    for f in files:
        all_hits.append(pd.read_csv(f))

    return all_hits

@pytest.fixture(scope="module")
def local_blast_output():
    """
    These csv files are output from topiary.ncbi._blast_blast._thread_manager
    """

    dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.abspath(os.path.join(dir,"data","local-blast-output"))

    files = glob.glob(os.path.join(base_dir,"*.csv"))
    files.sort()

    all_hits = []
    for f in files:
        all_hits.append(pd.read_csv(f))

    return all_hits

@pytest.fixture(scope="module")
def seed_dataframes():

    dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.abspath(os.path.join(dir,"data","seed-dataframes"))
    files = os.listdir(base_dir)

    out_dict = {}
    for f in files:
        out_file = os.path.join(base_dir,f)
        if os.path.isfile(out_file):
            out_dict[f] = out_file

    return out_dict

@pytest.fixture(scope="module")
def user_seed_dataframes():

    dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.abspath(os.path.join(dir,
                                            "data",
                                            "seed-dataframes",
                                            "user-seed-dataframes"))


    files = glob.glob(os.path.join(base_dir,"*.xlsx"))

    out_dict = {}
    for f in files:
        out_dict[os.path.split(f)[-1]] = f

    return out_dict



@pytest.fixture(scope="module")
def esummary_assembly_records():

    dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.abspath(os.path.join(dir,"data","ncbi-assembly-esummary-output"))
    json_files = glob.glob(os.path.join(base_dir,"*.json"))

    out_dict = {}
    for json_file in json_files:
        key = os.path.basename(json_file).split(".")[0]
        f = open(json_file,'r')
        out_dict[key] = json.load(f)
        f.close()

    return out_dict

@pytest.fixture(scope="module")
def make_blast_db_files():

    dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.abspath(os.path.join(dir,"data","make_blast_db_files"))
    files = glob.glob(os.path.join(base_dir,"*"))

    out_dict = {}
    for f in files:
        key = os.path.basename(f)
        out_dict[key] = f

    return out_dict

@pytest.fixture(scope="module")
def for_real_inference():

    dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.abspath(os.path.join(dir,"data","for-real-inference"))
    files = glob.glob(os.path.join(base_dir,"*"))

    out_dict = {}
    for f in files:
        key = os.path.basename(f)
        out_dict[key] = f

    return out_dict

@pytest.fixture(scope="module")
def df_with_species_not_resolvable():

    dir = os.path.dirname(os.path.realpath(__file__))
    df_file = os.path.abspath(os.path.join(dir,"data","test_dataframes","df-with-species-not-resolvable.csv"))
    df = pd.read_csv(df_file)

    return df

@pytest.fixture(scope="module")
def ftp_test_files():

    dir = os.path.dirname(os.path.realpath(__file__))

    base_dir = os.path.abspath(os.path.join(dir,"data","ftp"))
    files = glob.glob(os.path.join(base_dir,"*"))

    out_dict = {}
    for f in files:
        key = os.path.basename(f)
        out_dict[key] = f

    return out_dict

@pytest.fixture(scope="module")
def generax_data():

    dir = os.path.dirname(os.path.realpath(__file__))

    base_dir = os.path.abspath(os.path.join(dir,"data","generax"))
    files = glob.glob(os.path.join(base_dir,"*"))

    out_dict = {}
    for f in files:
        key = os.path.basename(f)
        out_dict[key] = f

    return out_dict

@pytest.fixture(scope="module")
def tiny_phylo():
    return get_files(os.path.join("data","tiny-phylo"))

@pytest.fixture(scope="module")
def small_phylo():
    return get_files(os.path.join("data","small-phylo"))

@pytest.fixture(scope="module")
def raxml_output():
    return get_files(os.path.join("data","raxml-output"))

@pytest.fixture(scope="module")
def alignment_input():
    return get_files(os.path.join("data","alignment-input"))
