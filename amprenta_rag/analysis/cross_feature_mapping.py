from typing import Dict, List
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.models import feature_pathway_map, gene_protein_map


def map_gene_to_proteins(gene_symbol: str, db: Session) -> List[str]:
    query = db.execute(gene_protein_map.select().where(gene_protein_map.c.gene_symbol == gene_symbol))
    return [row.uniprot_id for row in query.fetchall()]


def map_protein_to_genes(uniprot_id: str, db: Session) -> List[str]:
    query = db.execute(gene_protein_map.select().where(gene_protein_map.c.uniprot_id == uniprot_id))
    return [row.gene_symbol for row in query.fetchall()]


def get_cross_omics_feature_neighbors(feature_id: UUID, db: Session) -> Dict:
    # Returns dict of related genes, proteins, pathways for a feature_id
    neighbors = {}
    gene_rows = db.execute(feature_pathway_map.select().where(feature_pathway_map.c.feature_id == feature_id))
    neighbors["pathways"] = [row.pathway_id for row in gene_rows.fetchall()]
    # Add other mappings as your schema enables (e.g., via joined query to gene_protein_map)
    return neighbors
