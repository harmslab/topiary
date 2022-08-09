
import pytest

import topiary
from topiary.ncbi.entrez.proteome import _get_genome_url
from topiary.ncbi.entrez.proteome import get_proteome_ids
from topiary.ncbi.entrez.proteome import get_proteome

import numpy as np
import pandas as pd

import datetime, os

def test__get_genome_url(esummary_assembly_records):

    # If parsing and sorting is doing what we think, we should get the following
    # out of our stack of urls for each json
    expected_output = {"Yersinia_pestis":"GCF_000222975.1_ASM22297v1",
                        "Pseudotsuga_menziesii":"GCA_001517045.1_DougFir1.0",
                        "Saccharomyces_cerevisiae":"GCA_000976875.2_Sc_YJM1244_v1",
                        "Schizosaccharomyces_pombe":"GCF_000002945.1_ASM294v2",
                        "Escherichia_coli":"GCF_023598925.1_ASM2359892v1",
                        "Salmonella_enterica":"GCF_023572745.1_ASM2357274v1",
                        "Thermus_thermophilus":"GCF_000091545.1_ASM9154v1",
                        "Homo_sapiens":"GCF_000001405.40_GRCh38.p14",
                        "Mus_musculus":"GCF_000001635.27_GRCm39",
                        "Gallus_gallus":"GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b",
                        "Danio_rerio":"GCF_000002035.6_GRCz11",
                        "Monosiga_brevicollis":"GCF_000002865.3_V1.0",
                        "Arabidopsis_thaliana":"GCF_000001735.4_TAIR10.1",
                        "Zea_mays":"GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0"}

    for k in esummary_assembly_records:
        output = []
        for r in esummary_assembly_records[k]:
            out = _get_genome_url(r)
            assert out is None or len(out) == 4
            if out:
                assert type(out[2]) is datetime.datetime
                output.append(out)

        output.sort()
        acc = output[0][3].split("/")[-1]

        assert acc == expected_output[k]

def test_get_proteome_ids():

    with pytest.raises(ValueError):
        get_proteome_ids(species=None,taxid=None)

    with pytest.raises(ValueError):
        get_proteome_ids(species="Homo sapiens",taxid="9606")

    # Bad taxid query
    id_list, err = get_proteome_ids(taxid=999999999999)
    assert id_list is None
    assert isinstance(err,str)

    # Bad species query
    id_list, err = get_proteome_ids(species="Not a species")
    assert id_list is None
    assert isinstance(err,str)

    # Human, by taxid
    by_taxid, err = get_proteome_ids(taxid=9606)
    assert isinstance(by_taxid,list)
    assert len(by_taxid) > 0
    assert err is None

    # Human, by species
    by_species, err = get_proteome_ids(species="Homo sapiens")

    # Make sure species and taxid queries bring down equivalent ids
    assert np.array_equal(by_species,by_taxid)


def test_get_proteome(tmpdir):

    cwd = os.getcwd()
    os.chdir(tmpdir)

    with pytest.raises(ValueError):
        get_proteome(species=None,taxid=None)

    with pytest.raises(ValueError):
        get_proteome(species="Homo sapiens",taxid="9606")

    # Actually pull down the human proteome, making sure it comes down and is
    # written to output.
    output1 = get_proteome(taxid=9606)
    assert os.path.isfile(output1)

    os.chdir(cwd)
