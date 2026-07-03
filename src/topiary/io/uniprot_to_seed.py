"""
Function for converting a UniProt FASTA file into a topiary seed dataframe.
"""

import topiary
from topiary.io import read_seed

import pandas as pd
import re
import io

def uniprot_to_seed(fasta_file,
                    read_with_read_seed=True):
    """
    Convert a UniProt FASTA file into a topiary seed dataframe.

    Parameters
    ----------
    fasta_file : str
        Path to the UniProt FASTA file.
    read_with_read_seed : bool, default=True
        Whether or not to pass the resulting dataframe through
        topiary.io.read_seed to ensure it is valid.

    Returns
    -------
    seed_df : pandas.DataFrame
        A topiary seed dataframe.
    """

    # Read the fasta file
    with open(fasta_file, "r") as f:
        content = f.read()

    # Split into entries
    entries = content.split(">")[1:]

    # To ensure global alias uniqueness (case-insensitive) across the whole
    # seed dataframe (to avoid ValueError in read_seed)
    seen_aliases = set()

    data = []
    for entry in entries:
        lines = entry.split("\n")
        header = lines[0]
        sequence = "".join([line.strip() for line in lines[1:]])

        # Parse UniProt header
        # >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
        
        # Split by pipe for db, accession, entry name
        parts = header.split("|", 2)
        if len(parts) < 3:
            # Not a standard UniProt header, skip or handle differently?
            # For now, skip
            continue
        
        db = parts[0]
        accession = parts[1]
        rest = parts[2]
        
        # Entry name is the first word of the rest
        entry_parts = rest.split(None, 1)
        entry_name = entry_parts[0]
        remainder = entry_parts[1] if len(entry_parts) > 1 else ""
        
        # Parse KV pairs: OS=, OX=, GN=, PE=, SV=
        kv_pattern = re.compile(r"(?:\s+|^)(OS|OX|GN|PE|SV)=")
        kv_matches = list(kv_pattern.finditer(remainder))
        
        if not kv_matches:
            protein_name = remainder.strip()
            kv_data = {}
        else:
            protein_name = remainder[:kv_matches[0].start()].strip()
            kv_data = {}
            for i, match in enumerate(kv_matches):
                key = match.group(1)
                start = match.end()
                if i + 1 < len(kv_matches):
                    end = kv_matches[i+1].start()
                else:
                    end = len(remainder)
                value = remainder[start:end].strip()
                kv_data[key] = value

        # Map to topiary columns
        species = kv_data.get("OS", "unknown")
        
        # Name is GN if available, else entry_name
        name = kv_data.get("GN", entry_name)

        # Track the name itself as a seen alias to prevent other proteins
        # from claiming it as an alias.
        seen_aliases.add(name.lower().strip())
        
        # Aliases - deduplicate globally (case-insensitive). Also ensure
        # that aliases do not collide with the name.
        raw_aliases = [accession, entry_name, protein_name]
        unique_for_this_protein = []
        for a in raw_aliases:
            if not a:
                continue
            
            a_lower = a.lower().strip()
            if a_lower not in seen_aliases:
                unique_for_this_protein.append(a)
                seen_aliases.add(a_lower)
        
        aliases = ";".join(unique_for_this_protein)
        
        data.append({
            "species": species,
            "name": name,
            "aliases": aliases,
            "sequence": sequence
        })

    df = pd.DataFrame(data)

    if read_with_read_seed:
        # Use read_seed to ensure it's valid. This returns (df, key_species,
        # paralog_patterns, species_aware). We only want the df.
        df, _, _, _ = read_seed(df)

    return df
