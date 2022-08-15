
import topiary
from topiary._private.interface import gen_seed
from topiary._private import check

import os
import json
import time
import copy
import shutil
import glob
import datetime

class Supervisor:
    """
    Class that keeps track of job status, creates and reads stereotyped
    calculation directories, and writes out job metadata in a standard json
    format.

    A supervisor-created directory will have the following structure:

    .. code_block::

        calc_dir
            input/
            working/
            output/
            run_parameters.json

    A Supervisor's primary job is to wrap, update, and allow programmatic access
    to the calculation attributes stored in run_parameters.json. This dictionary
    has the following entries. Only the "calc_status", "version", and "seed"
    keys are gauranteed to exist.

    + "calc_status" ("empty","running","crashed","complete")
    + "version" (topiary version)
    + "seed" (random integer seed)
    + "creation_time" (when directory was made via create_calc_dir)
    + "completion_time" (when finalize was called)
    + "calc_type" (what type of calculation this is)
    + "model" (phylogenetic model)
    + "events" (list of events passed to supervisor during run)
      The "events" entries will each be dictionaries with the following sorts of
      keys:

      + "description" (description of what the calculation is)
      + "local_directory" (directory in which job was run, relative to calc_dir)
      + "time" (time this was recorded)
      + **kwargs (other keys passed in as meta data)

    + "previous_entries" (This is a list holding run_parameters of all previous
      calculations managed by the supervisor. These dictionaries will have a
      "calc_dir" key added that keeps track of directory name they came from.
      The list is ordered earliest to latest).

    These parameters will are written out to run_parameters.json.

    This object also exposes some convenience properties. These are absolute
    paths to the corresponding file in the input directory. If the file does not
    exist, these properties are None.

    + gene_tree (input/gene-tree.newick)
    + species_tree (input/species-tree.newick)
    + reconciled_tree (input/reconciled_tree.newick)
    + alignment (input/alignment.py)

    Other important properties are:

    + df (topiary dataframe corresponding to input/dataframe.csv)
    + model (string model)
    + seed (seed)

    Note this is **not** thread-safe.
    """

    keep_on_increment = ["model","version","seed"]

    def __init__(self,calc_dir=None,seed=None):
        """
        Create a Supervisor instance.

        Parameters
        ----------
        calc_dir : str, optional
            load the supervisor from an existing directory
        seed : int, optional
            seed for calcultions. if not specified, generate a random seed
        """

        # seed argument is overloaded. Interpret based on type
        if seed is None:
            seed = gen_seed()

        # By this point, seed better be an int ...
        seed = check.check_int(seed,minimum_allowed=0)

        self._run_parameters = {"calc_status":"empty",
                                "version":topiary.__version__,
                                "seed":seed}
        self._calc_dir = calc_dir
        self._starting_dir = os.path.abspath(os.getcwd())
        self._df = None
        self._run_parameters["events"] = []

        if self._calc_dir is not None:
            self._calc_dir = os.path.abspath(self._calc_dir)
            if not os.path.exists(self._calc_dir):
                err = f"\nThe directory '{calc_dir}' does not exist. calc_dir\n"
                err += "should be a previously-run calculation. To create a\n"
                err += "calculation directory, please use sv = Supervisor()\n"
                err += f"then sv.create_calc_dir('{calc_dir}')\n\n"
                raise ValueError(err)

            self._load_existing()


    def _load_existing(self):
        """
        Load information about an existing calculation.
        """

        # Grab the run parameters
        try:
            f = open(os.path.join(self._calc_dir,"run_parameters.json"),'r')
            run_parameters = json.load(f)
            f.close()

        except FileNotFoundError as e:
            err = f"\nCould not read previous directory '{self._calc_dir}'. This\n"
            err += "directory should contain run_parameters.json. This file was\n"
            err += "not found.\n\n"
            raise FileNotFoundError(err) from e

        except json.JSONDecodeError as e:

            err = f"\nCould not read previous directory '{self._calc_dir}'. This\n"
            err += "directory should contain run_parameters.json. This file exists\n"
            err += "but could not be read!\n\n"
            raise ValueError(err) from e

        self._run_parameters = run_parameters

        # Read df_file if it exists
        df_file = os.path.join(self._calc_dir,"input","dataframe.csv")
        if os.path.isfile(df_file):
            self._df = topiary.read_dataframe(df_file)


    def _increment(self,new_calc_dir):
        """
        Increment the calculation with creation of a new directory, moving
        current calc into "previous_entries" list and popping out most
        attributes so we are reset for the next calculation.

        Parameters
        ----------
        new_calc_dir : str
            new calculation directory
        """

        if self._run_parameters["calc_status"] == "empty":
            self._calc_dir = os.path.abspath(new_calc_dir)
            return

        keys = [k for k in self._run_parameters.keys() if k != "previous_entries"]
        if "previous_entries" not in self._run_parameters:
            self._run_parameters["previous_entries"] = []

        # Create previous entry
        self._run_parameters["previous_entries"].append({})
        prev = self._run_parameters["previous_entries"][-1]

        # Record calc and start directories
        prev["calc_dir"] = self._calc_dir
        prev["starting_dir"] = self._starting_dir

        for k in keys:
            prev[k] = copy.deepcopy(self._run_parameters[k])
            if k not in self.keep_on_increment:
                self._run_parameters.pop(k)

        self._run_parameters["events"] = []

        self._calc_dir = os.path.abspath(new_calc_dir)

    def create_calc_dir(self,
                        calc_dir,
                        calc_type,
                        overwrite=False,
                        force=False,
                        df=None,
                        gene_tree=None,
                        species_tree=None,
                        reconciled_tree=None,
                        model=None):
        """
        Create a calculation directory in a stereotyped way.

        Parameters
        ----------
        calc_dir : str
            directory to create
        calc_type : str
            type of calculation
        overwrite : bool, default=False
            whether or not to overwrite existing calc_dir
        force : bool, default=False
            force creating new calculation even if last calculation did not
            complete (will likely crash)
        df : pandas.DataFrame or str, optional
            dataframe to use for calculation (goes into input/dataframe.csv).
            Overwrites whatever came from previous directory.
        gene_tree : str, ete3.Tree, dendropy.tree, optional
            gene_tree file for calculation (goes into input/gene-tree.newick).
            If this an ete3 or dendropy tree, it will be written out with leaf
            names and branch lengths; all other data will be dropped.
        species_tree : str, ete3.Tree, dendropy.tree, optional
            species_tree file for calculation (goes into input/species-tree.newick).
            If this an ete3 or dendropy tree, it will be written out with leaf
            names; all other data will be dropped.
        reconciled_tree : str, ete3.Tree, dendropy.tree, optional
            species_tree file for calculation (goes into input/reconciled-tree.newick).
            If this an ete3 or dendropy tree, it will be written out with leaf
            names; all other data will be dropped.
        model : str, optional
            phylogenetic model recognized by raxml and generax
        """

        if self.status not in ["complete","empty"]:
            if not force:
                err = f"previous calculation has status '{self.status}', not\n"
                err += f"'complete'. Cannot create '{calc_dir}'. To proceed anyway,\n"
                err += "set force = True.\n\n"
                raise ValueError(err)

        # Make sure we can make the calc_dir
        if os.path.exists(calc_dir):
            if overwrite:
                shutil.rmtree(calc_dir)
            else:
                err = f"\ncalc_dir '{calc_dir}' already exists\n\n"
                raise FileExistsError(err)

        self._starting_dir = os.path.abspath(os.getcwd())

        # Create new directory
        os.mkdir(calc_dir)
        os.mkdir(os.path.join(calc_dir,"input"))
        os.mkdir(os.path.join(calc_dir,"working"))
        os.mkdir(os.path.join(calc_dir,"output"))

        # Increment internal run_parameters tracker. (Note this will set
        # self._calc_dir)
        self._increment(calc_dir)

        # Get ready to load in previous data
        self._df = None
        self._run_parameters["gene_tree"] = None
        self._run_parameters["species_tree"] = None
        self._run_parameters["reconciled_tree"] = None
        self._run_parameters["alignment"] = None
        input_dir = os.path.join(self._calc_dir,"input")

        # ---------------------------------------------------------------------
        # Load df and tree from previous run

        # See if there is a previous output directory to read
        if "previous_entries" in self._run_parameters:

            last = self._run_parameters["previous_entries"][-1]
            previous_dir = os.path.join(last["calc_dir"],"output")

            prev_df = os.path.join(previous_dir,"dataframe.csv")
            if os.path.isfile(prev_df):
                input_df = os.path.join(input_dir,"dataframe.csv")
                shutil.copy(prev_df,input_df)
                self._df = topiary.read_dataframe(input_df)

            trees = ["gene-tree.newick","species-tree.newick","reconciled-tree.newick"]
            for t in trees:
                prev_tree = os.path.join(previous_dir,t)
                if os.path.isfile(prev_tree):
                    input_tree = os.path.join(input_dir,t)
                    shutil.copy(prev_tree,input_tree)
                    self._run_parameters["_".join(t.split(".")[0].split("-"))] = input_tree

        # ---------------------------------------------------------------------
        # Load arguments (after we've already processed stuff from previous
        # directory so we overwrite with arguments)

        if df is not None:
            df = topiary.read_dataframe(df)
            topiary.write_dataframe(df,os.path.join(input_dir,"dataframe.csv"),
                                    overwrite=True)
            self._df = df

        tree_names = ["gene-tree","species-tree","reconciled-tree"]
        for i, tree in enumerate([gene_tree,species_tree,reconciled_tree]):

            if tree is not None:

                # Make sure it is readable (even if a string/file)
                T = topiary.io.read_tree(tree)

                out_tree = os.path.join(input_dir,f"{tree_names[i]}.newick")
                if os.path.isfile(str(tree)):
                    shutil.copy(tree,out_tree)
                else:
                    T.write(outfile=out_tree,format=5)

                self._run_parameters["_".join(tree_names[i].split("-"))] = out_tree

        if model is not None:
            self._run_parameters["model"] = model

        # ---------------------------------------------------------------------
        # Write alignment.phy from dataframe

        if self._df is not None:
            aln_file = os.path.join(input_dir,"alignment.phy")
            topiary.io.write_phy(self._df,aln_file)
            self._run_parameters["alignment"] = aln_file

        # ---------------------------------------------------------------------
        # Record calc status

        self._run_parameters["calc_type"] = calc_type
        self._run_parameters["calc_status"] = "running"
        self._run_parameters["creation_time"] = time.time()
        self.write_json()

        # ---------------------------------------------------------------------
        # Copy input files we are going to want in output

        # First, get newick files from previous output directory and put in new
        # output. These might get wiped out by new files, but we want to make
        # sure everything comes in appropriately.
        self.copy_output_to_output("*.newick")

        # Next get newick files from the input directory.
        for n in glob.glob(os.path.join(self.input_dir,"*.newick")):
            self.stash(n)

        # Finaly get input/dataframe.csv if one was passed in
        try:
            self.stash(os.path.join(self.input_dir,"dataframe.csv"))
        except FileNotFoundError:
            pass


        # ---------------------------------------------------------------------
        # Print status

        d = calc_dir
        if not d[0] == "/":
            d = f"./{d}"
        out = ["\n----------------------------------------------------------------------\n"]
        out.append(f"topiary is starting a {self.calc_type} calculation in {d}:\n")
        print("\n".join(out),flush=True)



    def check_required(self,
                      required_files=[],
                      required_values=[]):
        """
        Check for required files and values for a particular calculation.

        Parameters
        ----------
        required_files : list, default=[]
            list of files to look for (basename; will look in input)
        required_values : list, default=[]
            list of values to look for in self._run_parameters. considered
            missing if value is None. Does not check range/type, just presence
            or None.
        """

        # Check for missign files
        missing = []
        for f in required_files:
            look_for = os.path.join(self.input_dir,f)
            if not os.path.isfile(look_for):
                missing.append(f)

        if len(missing) > 0:
            err = "\nSome necessary files are missing:\n\n"
            for m in missing:
                err += f"    {m}\n"
            err += ""
            raise FileNotFoundError(err)

        # Check for missing values
        missing = []
        for r in required_values:
            try:
                value = self._run_parameters[r]
                if value is None:
                    raise KeyError
            except KeyError:
                missing.append(r)

        if len(missing) > 0:
            err = "\nSome required values are not defined:\n\n"
            for m in missing:
                err += f"    {m}\n"
            err += ""
            raise ValueError(err)

    def copy_output_to_output(self,pattern="*.newick"):
        """
        Copy files from the previous output directory and put them into the new
        output directory. This should typically be run *before* any newly
        calculated outputs are stashed, so the new files replace the old ones.

        Parameters
        ----------
        pattern : str
            glob pattern for the file(s) to copy
        """

        if self.previous_entries is not None:
            prev_dir = self.previous_entries[-1]["calc_dir"]
            existing_files = glob.glob(os.path.join(prev_dir,"output",pattern))
            for f in existing_files:
                filename = os.path.basename(f)
                shutil.copy(f,os.path.join(self.output_dir,filename))

    def stash(self,to_stash,target_name=None,target_dir="output"):
        """
        Copy the file or directory to_stash into the run directory. This will
        overwrite without warning.

        Parameters
        ----------
        to_stash : str
            file or directory to copy in
        target_name : str, optional
            what to name file or directory that is being stored. If not specified,
            use basename. If this target has directories, these will be created.
        target_dir : str, default="output"
            directory in topiary in which to stash (should be "input","working",
            or "output")

        Returns
        -------
        final_target_path : str
            absolute path to the newly stashed file or directory
        """

        if not os.path.exists(to_stash):
            err = f"\nto_stash '{to_stash}' does not exist.\n\n"
            raise FileNotFoundError(err)

        if target_dir not in ["input","working","output"]:
            err = f"target_dir is '{target_dir}'. Should be 'input', 'working'\n"
            err += "or 'output'.\n\n"
            raise ValueError(err)

        # Make target_dir absolute
        target_dir = os.path.join(self._calc_dir,target_dir)

        # Get target_name
        if target_name is None:
            target_name = os.path.basename(to_stash)

        # Append directories, if necessary, to wherever we are stashing this
        target_split = os.path.split(target_name)
        for t in target_split[:-1]:
            t_base = os.path.join(target_dir,t)

            # Nuke files that share path in directory
            if os.path.exists(t_base) and not os.path.isdir(t_base):
                os.remove(t_base)

            if not os.path.exists(t_base):
                os.mkdir(t_base)

        # Final path
        final_target_path = os.path.join(target_dir,target_name)

        if os.path.isdir(to_stash):
            shutil.copytree(to_stash,final_target_path)
        else:
            shutil.copy(to_stash,final_target_path)

        return final_target_path


    def write_json(self):
        """
        Write out a json with the current status.
        """

        f = open(os.path.join(self._calc_dir,"run_parameters.json"),"w")
        json.dump(self._run_parameters,f)
        f.close()


    def update(self,key,value):
        """
        Update run_parameters with a key/value pair. Used to store parameter
        values and meta data (cutoffs, user options, etc.) for each calculation.
        Will not function after self.finalize() has been called.

        Parameters
        ----------
        key : str
            parameter to update
        value :
            value to store
        """

        if self.status in ["complete","crashed"]:
            err = "\nCannot set values after the calculation is finalized.\n\n"
            raise RuntimeError(err)

        if not issubclass(type(key),str):
            err = "\njson keys must be strings.\n\n"
            raise ValueError(err)

        self._run_parameters[key] = value

    def event(self,description,**kwargs):
        """
        Start an individual calculation.

        Parameters
        ----------
        local_directory : str
            directory (relative to starting_dir) in which a calculation is
            running. This is relative if the directory is in a supervisor
            directory. It will be absolute if the current working directory
            is outside.
        description : str, optional
            string description of what is being gone
        **kwargs : kwargs useful for keeping track of whatever is being run
        """

        if self._run_parameters["calc_status"] != "running":
            s = self._run_parameters["calc_status"]
            err = f"supervisor status is '{s}'. Status must be 'running' to\n"
            err += "record an event.\n\n"
            raise RuntimeError(err)

        try:
            self._run_parameters["events"]
        except KeyError:
            self._run_parameters["events"] = []

        # Get path relative to calc_dir where this calculation is running
        current_dir = os.getcwd()
        os.chdir(self.calc_dir)
        local_directory = os.path.relpath(current_dir)
        os.chdir(current_dir)

        # If outside of calc_dir, record as an absolute path
        if local_directory.startswith(".."):
            local_directory = current_dir

        this_time = time.time()
        self._run_parameters["events"].append({"description":description,
                                               "local_directory":local_directory,
                                               "time":this_time})

        for k in kwargs:
            if k in ["local_directory","time"]:
                err = "\nevent keys 'local_directory' and 'time' are reserved\n\n"
                raise ValueError(err)
            self._run_parameters["events"][-1][k] = kwargs[k]


        # Print output
        dt = this_time - self._run_parameters["creation_time"]
        pretty_time = f"{str(datetime.timedelta(seconds=dt))} (H:M:S)"
        print(f"{description}, {pretty_time}",flush=True)

        self.write_json()

    def finalize(self,successful=True,plot_if_success=False):
        """
        Finalize a total set of calculations. Once invoked .event() and .update()
        can no longer be run.

        Parameters
        ----------
        successful : bool
            whether or not job actually completed successfully

        Returns
        -------
        None OR typlot.canvas if plot requested
        """

        # Do not update if already finalized or empty
        if self.status in ["complete","crashed","empty"]:
            return None

        # Update json
        self._run_parameters["completion_time"] = time.time()
        if successful:
            self._run_parameters["calc_status"] = "complete"
        else:
            self._run_parameters["calc_status"] = "crashed"

            # Return to starting directory, whatever that was
            os.chdir(self.starting_dir)

        self.write_json()

        # Get time to run
        dt = self._run_parameters["completion_time"] - self._run_parameters["creation_time"]
        pretty_time = f"{str(datetime.timedelta(seconds=dt))} (H:M:S)"

        # Get pretty names for directories
        calc_dir = os.path.relpath(self.calc_dir)
        if calc_dir[0] != "/":
            calc_dir = os.path.join("./",calc_dir)

        output_dir = os.path.relpath(self.output_dir)
        if output_dir[0] != "/":
            output_dir = os.path.join("./",output_dir)

        working_dir = os.path.relpath(self.working_dir)
        if working_dir[0] != "/":
            working_dir = os.path.join("./",working_dir)

        out = ["\n"]
        out.append(f"topiary ran a {self.calc_type} calculation in {calc_dir}:\n")
        if successful:
            out.append(f"+ Completed in {pretty_time}")
            out.append(f"+ Wrote results to {output_dir}")
        else:
            out.append(f"+ Crashed after {pretty_time}")
            out.append(f"+ Please check {working_dir}")

        out.append("\n----------------------------------------------------------------------\n\n")

        print("\n".join(out),flush=True)

        # If plot requested, generate plot and return
        if successful and plot_if_success:
            return topiary.draw.tree(calculation=self.calc_dir,
                                     output_file=os.path.join(self.output_dir,
                                                              "summary-tree.pdf"))

        return None


    @property
    def status(self):
        """
        Current status of the calculation. Will be one of:

        + "empty" (no directory created or files loaded)
        + "running" (directory created and presumably some calculation running)
        + "crashed" (a calculation crashed)
        + "complete" (calculation was finalized successfully).
        """
        return self._run_parameters["calc_status"]


    @property
    def starting_dir(self):
        """
        absolute path to directory in which Supervisor() was instantiated or
        last time create_calc_dir() was invoked.
        """
        return self._starting_dir

    @property
    def calc_dir(self):
        """
        absolute path to calc_dir (directory containing) input/, working/,
        output/, run_parameters.json
        """

        if self._calc_dir is not None:
            return self._calc_dir
        else:
            return None

    @property
    def input_dir(self):
        """
        absolute path to calc_dir/input
        """
        if self._calc_dir is not None:
            return os.path.join(self._calc_dir,"input")
        else:
            return None

    @property
    def working_dir(self):
        """
        absolute path to calc_dir/working
        """
        if self._calc_dir is not None:
            return os.path.join(self._calc_dir,"working")
        else:
            return None

    @property
    def output_dir(self):
        """
        absolute path to calc_dir/output
        """
        if self._calc_dir is not None:
            return os.path.join(self._calc_dir,"output")
        else:
            return None

    @property
    def calc_type(self):
        """
        Type of calculation (as reported by run_parameters["calc_type"])
        """
        if "calc_type" in self._run_parameters:
            return self._run_parameters["calc_type"]
        else:
            return None

    @property
    def gene_tree(self):
        """
        Absolute path to input gene tree for calculation (as reported by
        run_parameters["gene_tree"])
        """

        if "gene_tree" in self._run_parameters:
            return self._run_parameters["gene_tree"]
        else:
            return None


    @property
    def species_tree(self):
        """
        Absolute path to input species tree for calculation (as reported by
        run_parameters["species_tree"])
        """

        if "gene_tree" in self._run_parameters:
            return self._run_parameters["species_tree"]
        else:
            return None

    @property
    def reconciled_tree(self):
        """
        Absolute path to input species tree for calculation (as reported by
        run_parameters["reconciled_tree"])
        """

        if "reconciled_tree" in self._run_parameters:
            return self._run_parameters["reconciled_tree"]
        else:
            return None


    @property
    def alignment(self):
        """
        Absolute path to input alignment file for calculation (as reported by
        run_parameters["alignment"])
        """
        if "alignment" in self._run_parameters:
            return self._run_parameters["alignment"]
        else:
            return None

    @property
    def seed(self):
        """
        Calculation seed (as reported by run_parameters["seed"])
        """
        return self._run_parameters["seed"]

    @property
    def model(self):
        """
        Calculation model (as reported by run_parameters["model"])
        """
        if "model" in self._run_parameters:
            return self._run_parameters["model"]
        else:
            return None

    @property
    def df(self):
        """
        Topiary instance corresponding to self.input_dir/dataframe.csv.
        """

        return self._df

    @property
    def previous_entries(self):
        """
        List of run_parameters dictionaries for all previous runs done in this
        calculation stack. The last entry is the most recent.
        """

        if "previous_entries" in self._run_parameters:
            return self._run_parameters["previous_entries"]
        else:
            return None

    @property
    def run_parameters(self):
        """
        run_parameters dictionary (equivalent to run_parameters.json).
        """

        return self._run_parameters

    @property
    def tree_prefix(self):
        """
        Class of tree from this output. Will be "reconciled", "gene",
        or None.
        """

        try:
            calc_type = self._run_parameters["calc_type"]
        except KeyError:
            return None

        if calc_type.split("_")[0] == "ml":
            return "gene"

        if calc_type.split("_")[0] == "reconcile":
            return "reconciled"

        return None
