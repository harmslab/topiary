
import topiary

from topiary.external.opentree import phylo_to_taxid

def test_init():

    # Make sure the phylo_to_taxid database is being read properly with a few
    # spot checks. 
    assert phylo_to_taxid["Forams"] == 29178
    assert phylo_to_taxid["Excavata"] == (2611352,2611341,136087)
    assert phylo_to_taxid["Life"] is None
    assert phylo_to_taxid["LIFE"] is None
    assert phylo_to_taxid["life"] is None
    assert phylo_to_taxid["Cnidarians"] == 6073
