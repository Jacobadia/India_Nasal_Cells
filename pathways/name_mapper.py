class GeneNameMapper:

    default_mapping_file = "../artifacts/genetype_lookup.txt"

    def __init__(self, mapping_file: str = default_mapping_file):
        self.mapping = self.load_mapping(mapping_file)

    def load_mapping(self, mapping_file: str) -> dict[str, str]:
        mapping: dict[str, str] = {}
        with open(mapping_file, "rt", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                assert len(parts) >= 2
                gene_id, gene_name = parts[0], parts[1]
                gene_id = gene_id.strip()
                gene_name = gene_name.strip()
                mapping[gene_id] = gene_name
        return mapping
    
    def map_gene_id(self, gene_id: str) -> str:
        if gene_id not in self.mapping:
            raise ValueError(f"Gene ID not found in mapping: {gene_id}")
        gene_name = self.mapping[gene_id]
        assert gene_name
        if gene_name.startswith("ENSG"):
            print(f"INFO: Gene ID {gene_id} maps to another gene ID {gene_name}, not a name.")
        return gene_name