
import pytest

import topiary
from topiary.external.ncbi.entrez.proteome import get_proteome, _get_genome_url

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

def test_get_proteome(tmpdir):

    with pytest.raises(ValueError):
        get_proteome(species=None,taxid=None)

    with pytest.raises(ValueError):
        get_proteome(species="Homo sapiens",taxid="9606")

    with pytest.raises(ValueError):
        get_proteome(taxid="9606",output_dir="NOT_REALLY_A_DIR")

    # Actually pull down the human proteome, both using taxid and species.
    # Make sure these get the same file.
    output1 = get_proteome(taxid=9606,output_dir=tmpdir)
    assert os.path.isfile(output1)
    output2 = get_proteome(species="Homo sapiens",output_dir=tmpdir)
    assert output1 == output2
