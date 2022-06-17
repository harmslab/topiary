__description__ = \
"""
Use entrez to download a proteome off of the NCBI.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-06-09"

import topiary
from topiary import _arg_processors

from Bio import Entrez

import os, urllib, datetime, re

def _get_genome_url(record):
    """
    Get the url pointing to a genome with information for selecting the best
    genome.

    If a url is found, returns a tuple where entries are:

        RefSeq_category:
            0 -> reference genome
            1 -> representative genome
            2 -> something else
        Url class:
            0 -> FtpPath_RefSeq
            1 -> FtpPath_GenBank
        LastUpdateDate (as datetime.datetime inference)
        url

    If one sorts a list of these tuples, the first entry will be the best
    genome. If no genome is found, function returns None.

    Parameters
    ----------
        record: esummary record returned by Entrez.read

    Return
    ------
        genome info or None
    """

    raw_date = re.sub("/","-",record["LastUpdateDate"])
    last_date = datetime.datetime.fromisoformat(raw_date)

    if record["RefSeq_category"] == "reference genome":

        for i, r in enumerate(["FtpPath_RefSeq","FtpPath_GenBank"]):
            if record[r] != "":
                this_url = str(record[r])
                return (0,i,last_date,this_url)

    elif record["RefSeq_category"] == "representative genome":

        for i, r in enumerate(["FtpPath_RefSeq","FtpPath_GenBank"]):
            if record[r] != "":
                this_url = str(record[r])
                return (1,i,last_date,this_url)

    else:

        for i, r in enumerate(["FtpPath_RefSeq","FtpPath_GenBank"]):
            if record[r] != "":
                this_url = str(record[r])
                return (2,i,last_date,this_url)

    return None

def get_proteome(taxid=None,species=None,output_dir="."):
    """
    Download a proteome from the NCBI.

    Parameters
    ----------
        taxid: taxid (integer or string version of the integer). Incompatible
               with `species` argument.
        species: bionomial name (i.e. Mus musculus). Incompatible with `taxid`
                 argument.
        output_dir: where to write the file locally

    Return
    ------
        the file we downloaded or None if no file downloaded
    """

    if species is None and taxid is None:
        err = "\nYou must specify either species or taxid, but not both.\n\n"
        raise ValueError(err)

    if species is not None and taxid is not None:
        err = "\nYou must specify either species or taxid, but not both.\n\n"
        raise ValueError(err)

    output_dir = str(output_dir)
    if not os.path.isdir(output_dir) or not os.access(output_dir,os.W_OK):
        err = "\nlocal_path must exist and be writable.\n\n"
        raise ValueError(err)

    # Do query using taxid. This will throw an error if species is not sane,
    # so is a helpful check before querying ncbi.
    if species:
        print(f"Downloading proteome for species '{species}'",flush=True)
        taxid = topiary.ncbi.get_taxid(species)
    else:
        print(f"Downloading proteome for taxid '{taxid}'",flush=True)

    # Make sure the taxid is sane
    taxid = _arg_processors.process_int(taxid,"taxid")

    # Look for assembly ids for this taxid
    query_text = f"(txid{taxid}[ORGN])"
    esearch_handle = Entrez.esearch(db="assembly",
                                    retmax=1000,
                                    term=query_text,
                                    idtype="acc")
    search_record = Entrez.read(esearch_handle)

    # Get assembly ids
    try:
        returned_ids = list(search_record["IdList"])
    except KeyError:
        err = "\nThe Entrez.esearch query failed.\n\n"
        raise RuntimeError(err)

    # Make sure something came back
    if len(returned_ids) == 0:
        err = f"\nThe query '{query_text}' returned no assemblies.\n\n"
        raise RuntimeError(err)

    # Now get summary data for these records.
    esummary_query = ",".join(returned_ids)
    esummary_handle = Entrez.esummary(db="assembly",id=esummary_query)
    esummary_record = Entrez.read(esummary_handle)

    try:
        records = esummary_record["DocumentSummarySet"]["DocumentSummary"]
    except KeyError:
        err = "\nThe Entrez.esummary query failed.\n\n"
        raise RuntimeError(err)

    # Get list of urls from all records
    urls = []
    for record in records:
        ret = _get_genome_url(record)
        if ret is not None:
            urls.append(ret)
    urls.sort()

    # Go through the urls from highest to lowest quality and try to download
    # the proteome.
    success = False
    for u in urls:

        try:
            genome_url = u[3]
            refseq_name = os.path.basename(genome_url)
            out_file = f"{refseq_name}_protein.faa.gz"
            remote_file = f"{genome_url}/{out_file}"
            local_file = os.path.join(output_dir,out_file)
            urllib.request.urlretrieve(remote_file, local_file)
            success = True
            break

        except (urllib.error.URLError,urllib.error.HTTPError):
            continue

    # Try to delete random file that gets downloaded when we make this query.
    try:
        os.remove(os.path.join(output_dir,"esummary_assembly.dtd"))
    except FileNotFoundError:
        pass


    if success:
        return local_file

    return None
