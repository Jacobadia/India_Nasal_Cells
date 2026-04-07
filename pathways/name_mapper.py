from dataclasses import dataclass

@dataclass
class MappingEntry:
    gene_name: str
    gene_type: str

class GeneNameMapper:

    default_mapping_file = "../data/genetype_lookup.txt"

    def __init__(self, mapping_file: str = default_mapping_file):
        self.mapping: dict[str, MappingEntry] = self.load_mapping(mapping_file)

    def load_mapping(self, mapping_file: str) -> dict[str, MappingEntry]:
        mapping: dict[str, MappingEntry] = {}
        with open(mapping_file, "rt", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                assert len(parts) == 3
                gene_id, gene_name, gene_type = parts[0], parts[1], parts[2]
                gene_id = gene_id.strip()
                gene_name = gene_name.strip()
                gene_type = gene_type.strip()
                mapping[gene_id] = MappingEntry(gene_name, gene_type)
        return mapping
    
    def map_gene_id(self, gene_id: str) -> str:
        """Maps the ensemble id to the gene symbol"""
        if gene_id not in self.mapping:
            raise ValueError(f"Gene ID not found in mapping: {gene_id}")
        gene_name = self.mapping[gene_id].gene_name
        assert gene_name
        if gene_name.startswith("ENSG"):
            print(f"INFO: Gene ID {gene_id} maps to another gene ID {gene_name}, not a name.")
        return gene_name

    def map_gene_type(self, gene_id: str) -> str:
        """maps the ensemble id to the gene type"""
        if gene_id not in self.mapping:
            raise ValueError(f"Gene ID not found in mapping: {gene_id}")
        gene_type = self.mapping[gene_id].gene_type
        assert gene_type
        return gene_type