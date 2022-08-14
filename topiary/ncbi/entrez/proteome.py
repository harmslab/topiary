"""
Use entrez to download a proteome from the NCBI.
"""

import topiary
from topiary._private import check
from topiary.ncbi.entrez.download import ncbi_ftp_download


from Bio import Entrez

import os, urllib, datetime, re
import ftplib

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

def get_proteome_ids(taxid=None,species=None):
    """
    Query entrez to get a list of proteome ids that match a taxid or species.
    This will not raise an error on failure, but will instead return None with
    an error string.

    Parameters
    ----------
    taxid: int or str, optional
        NCBI taxid (integer or string version of the integer). Incompatible with
        `species` argument. At least taxid or species must be specified.
    species : str, optional
        bionomial name of species (i.e. Mus musculus). Incompatible with `taxid`
        argument. At least taxid or species must be specified.

    Returns
    -------
    returned_ids : list or None
        list of proteome ids. None if no ids found/error.
    err: str or None
        descriptive error if no returned_ids. None if no error.
    """

    if species is None and taxid is None:
        err = "\nYou must specify either species or taxid, but not both.\n\n"
        raise ValueError(err)

    if species is not None and taxid is not None:
        err = "\nYou must specify either species or taxid, but not both.\n\n"
        raise ValueError(err)

    # Do query using taxid. This will throw an error if species is not sane,
    # so is a helpful check before querying ncbi.
    if species:

        try:
            taxid = topiary.ncbi.get_taxid(species)
        except RuntimeError:
            err = f"\nCould not find taxid for species '{species}'\n\n"
            return None, err

    # Make sure the taxid is sane
    taxid = check.check_int(taxid,"taxid")

    returned_ids = None

    # Look for assembly ids for this taxid. Get reference genome first, then
    # go on to not reference. Sort
    filters = ['(latest[filter] AND "reference genome"[filter] AND all[filter] NOT anomalous[filter])',
               '(latest[filter] AND all[filter] NOT anomalous[filter])']
    for ref_filter in filters:

        query_text = f"txid{taxid}[ORGN] AND {ref_filter}"
        esearch_handle = Entrez.esearch(db="assembly",
                                        retmax=50,
                                        term=query_text,
                                        idtype="acc")
        search_record = Entrez.read(esearch_handle)

        # Get assembly ids
        try:
            returned_ids = list(search_record["IdList"])
        except KeyError:
            continue

        if len(returned_ids) == 0:
            continue

    if returned_ids is None:
        err = "\nThe Entrez.esearch query failed.\n\n"
        return None, err

    # Make sure something came back
    if len(returned_ids) == 0:
        err = f"\nThe Entrez.eserch query '{query_text}' returned no assemblies.\n\n"
        return None, err

    # Try to delete random file that gets downloaded when we make this query.
    try:
        os.remove("esummary_assembly.dtd")
    except FileNotFoundError:
        pass

    return returned_ids, None

def get_proteome(taxid=None,species=None):
    """
    Use entrez to download a proteome from the NCBI.

    Parameters
    ----------
    taxid: int or str, optional
        NCBI taxid (integer or string version of the integer). Incompatible with
        `species` argument. At least taxid or species must be specified.
    species : str, optional
        bionomial name of species (i.e. Mus musculus). Incompatible with `taxid`
        argument. At least taxid or species must be specified.

    Returns
    -------
    proteome_file : str or None
        the file we downloaded or None if no file downloaded
    """

    # Get proteome ids to download. Validate taxid and species arguments via
    # this function
    returned_ids, err = get_proteome_ids(taxid=taxid,species=species)
    if returned_ids is None:
        raise RuntimeError(err)

    if species is not None:
        print(f"Downloading proteome for species '{species}'",flush=True)
    else:
        print(f"Downloading proteome for taxid '{taxid}'",flush=True)


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

        genome_url = u[3]
        refseq_name = os.path.basename(genome_url)
        out_file = f"{refseq_name}_protein.faa.gz"
        remote_file = f"{genome_url}/{out_file}"

        try:
            ncbi_ftp_download(genome_url,file_base="_protein.faa.gz")
        except (ftplib.error_perm,RuntimeError):
            continue

        success = True
        break

    # Try to delete random file that can get downloaded when we make this query.
    try:
        os.remove("esummary_assembly.dtd")
    except FileNotFoundError:
        pass

    if success:
        return out_file

    return None
