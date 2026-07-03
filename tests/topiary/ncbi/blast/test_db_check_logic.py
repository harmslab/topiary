
import pytest
import os
from topiary.ncbi.blast.local import _prepare_for_blast
from topiary.ncbi.blast.make import make_blast_db

def test_db_check_logic(tmpdir):
    """
    Test the robust BLAST database existence check logic in _prepare_for_blast 
    and make_blast_db.
    """

    # Move into temporary directory
    os.chdir(tmpdir)

    # 1. Test single volume protein database (.psq)
    f = open("prot_single.psq","w")
    f.write("\n")
    f.close()

    # Should pass
    _prepare_for_blast(sequence="ACGT",
                       db="prot_single",
                       blast_program="blastp",
                       hitlist_size=10,
                       e_value_cutoff=0.001,
                       gapcosts=(11,1),
                       num_threads=1,
                       kwargs={},
                       test_skip_blast_program_check=True)

    # 2. Test multi-volume protein database with alias (.pal, .00.psq)
    f = open("prot_multi.pal","w")
    f.write("\n")
    f.close()
    f = open("prot_multi.00.psq","w")
    f.write("\n")
    f.close()

    # Should pass
    _prepare_for_blast(sequence="ACGT",
                       db="prot_multi",
                       blast_program="blastp",
                       hitlist_size=10,
                       e_value_cutoff=0.001,
                       gapcosts=(11,1),
                       num_threads=1,
                       kwargs={},
                       test_skip_blast_program_check=True)

    # 3. Test nucleotide database (.nsq)
    f = open("nuc_single.nsq","w")
    f.write("\n")
    f.close()

    # Should pass
    _prepare_for_blast(sequence="ACGT",
                       db="nuc_single",
                       blast_program="blastn",
                       hitlist_size=10,
                       e_value_cutoff=0.001,
                       gapcosts=(11,1),
                       num_threads=1,
                       kwargs={},
                       test_skip_blast_program_check=True)

    # 4. Test multi-volume nucleotide database with alias (.nal, .00.nsq)
    f = open("nuc_multi.nal","w")
    f.write("\n")
    f.close()
    f = open("nuc_multi.00.nsq","w")
    f.write("\n")
    f.close()

    # Should pass
    _prepare_for_blast(sequence="ACGT",
                       db="nuc_multi",
                       blast_program="blastn",
                       hitlist_size=10,
                       e_value_cutoff=0.001,
                       gapcosts=(11,1),
                       num_threads=1,
                       kwargs={},
                       test_skip_blast_program_check=True)

    # 5. Test non-existent database
    with pytest.raises(FileNotFoundError):
        _prepare_for_blast(sequence="ACGT",
                           db="non_existent",
                           blast_program="blastp",
                           hitlist_size=10,
                           e_value_cutoff=0.001,
                           gapcosts=(11,1),
                           num_threads=1,
                           kwargs={},
                           test_skip_blast_program_check=True)

def test_make_db_check_logic(tmpdir):
    """
    Test the robust BLAST database existence check logic in make_blast_db.
    """

    # Move into temporary directory
    os.chdir(tmpdir)

    # 1. Test existing .psq
    f = open("exists_psq.psq","w")
    f.write("\n")
    f.close()

    # Should raise FileExistsError because overwrite=False
    with pytest.raises(FileExistsError):
        make_blast_db(input_files=["test.faa"],db_name="exists_psq",overwrite=False)

    # 2. Test existing .pal
    f = open("exists_pal.pal","w")
    f.write("\n")
    f.close()

    # Should raise FileExistsError
    with pytest.raises(FileExistsError):
        make_blast_db(input_files=["test.faa"],db_name="exists_pal",overwrite=False)

    # 3. Test existing .00.psq
    f = open("exists_vol.00.psq","w")
    f.write("\n")
    f.close()

    # Should raise FileExistsError
    with pytest.raises(FileExistsError):
        make_blast_db(input_files=["test.faa"],db_name="exists_vol",overwrite=False)
