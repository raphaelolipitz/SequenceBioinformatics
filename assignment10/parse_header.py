from typing import Tuple, List

def parse_header(line: str) -> Tuple[str, str, str]:

    if (line[0] == '>' or line[0] == '+'):
        accession = line[1:line.index(' ')]
    else:
        accession = line[0:line.index(' ')]

    taxon = line[line.rindex('[') + 1:line.rindex(']')]

    protein = line[line.index(' ') + 1:line.index(' [')]

    return (accession, taxon, protein)