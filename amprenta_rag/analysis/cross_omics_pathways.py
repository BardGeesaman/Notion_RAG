from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Set
from uuid import UUID

from amprenta_rag.analysis.pathway.enrichment import perform_pathway_enrichment
from amprenta_rag.analysis.pathway.models import Pathway
from amprenta_rag.database.models import Dataset, Feature, Program, dataset_feature_assoc
from amprenta_rag.database.session import db_session


@dataclass
class OmicsFeatureSet:
    transcriptomics: List[str]
    proteomics: List[str]
    metabolomics: List[str]
    lipidomics: List[str]

    def all_features(self) -> Set[str]:
        return set(
            self.transcriptomics
            + self.proteomics
            + self.metabolomics
            + self.lipidomics
        )

    def feature_types(self) -> Set[str]:
        types = set()
        if self.transcriptomics:
            types.add("gene")
        if self.proteomics:
            types.add("protein")
        if self.metabolomics:
            types.add("metabolite")
        if self.lipidomics:
            types.add("lipid")
        return types

    def asdict(self) -> Dict[str, List[str]]:
        return asdict(self)


@dataclass
class ConvergentPathway:
    pathway_id: str
    name: str
    source: str
    adjusted_p_value: float
    enrichment_ratio: float
    matched_by_omics: Dict[str, List[str]]

    def asdict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class CrossOmicsEnrichmentResult:
    program_id: UUID
    pathways: List[ConvergentPathway]
    features: OmicsFeatureSet

    def asdict(self) -> Dict[str, object]:
        return {
            "program_id": self.program_id,
            "pathways": [p.asdict() for p in self.pathways],
            "features": self.features.asdict(),
        }


def _collect_features_for_datasets(dataset_ids: List[UUID]) -> OmicsFeatureSet:
    with db_session() as db:
        feats: List[Feature] = (
            db.query(Feature)
            .join(dataset_feature_assoc, Feature.id == dataset_feature_assoc.c.feature_id)
            .filter(dataset_feature_assoc.c.dataset_id.in_(dataset_ids))
            .all()
        )
    by_type: Dict[str, List[str]] = {
        "gene": [],
        "protein": [],
        "metabolite": [],
        "lipid": [],
    }
    for f in feats:
        ftype = getattr(f, "feature_type", None) or getattr(f, "type", None)
        name = getattr(f, "name", None)
        if not name:
            continue
        if ftype in by_type:
            by_type[ftype].append(name)
    return OmicsFeatureSet(
        transcriptomics=by_type["gene"],
        proteomics=by_type["protein"],
        metabolomics=by_type["metabolite"],
        lipidomics=by_type["lipid"],
    )


def analyze_pathway_convergence(features: OmicsFeatureSet) -> List[ConvergentPathway]:
    all_features = features.all_features()
    types = features.feature_types()
    if not all_features or not types:
        return []

    # Map feature -> omics list for matched_by_omics later
    feature_type_map: Dict[str, List[str]] = {}
    for f in features.transcriptomics:
        feature_type_map.setdefault(f, []).append("transcriptomics")
    for f in features.proteomics:
        feature_type_map.setdefault(f, []).append("proteomics")
    for f in features.metabolomics:
        feature_type_map.setdefault(f, []).append("metabolomics")
    for f in features.lipidomics:
        feature_type_map.setdefault(f, []).append("lipidomics")

    enrich_results = perform_pathway_enrichment(
        input_features=all_features,
        input_feature_types=types,
        p_value_threshold=0.05,
    )
    convergent: List[ConvergentPathway] = []
    for res in enrich_results:
        matched_by_omics: Dict[str, List[str]] = {
            "transcriptomics": [],
            "proteomics": [],
            "metabolomics": [],
            "lipidomics": [],
        }
        for m in res.matched_features:
            for ot in feature_type_map.get(m, []):
                matched_by_omics[ot].append(m)
        convergent.append(
            ConvergentPathway(
                pathway_id=res.pathway.pathway_id,
                name=res.pathway.name,
                source=res.pathway.source,
                adjusted_p_value=res.adjusted_p_value,
                enrichment_ratio=res.enrichment_ratio,
                matched_by_omics=matched_by_omics,
            )
        )
    return convergent


def combine_omics_features(dataset_ids: List[UUID]) -> OmicsFeatureSet:
    return _collect_features_for_datasets(dataset_ids)


def get_cross_omics_enrichment(program_id: UUID) -> CrossOmicsEnrichmentResult:
    with db_session() as db:
        program: Optional[Program] = db.query(Program).filter(Program.id == program_id).first()
        datasets = getattr(program, "datasets", []) if program else []
        dataset_ids = [d.id for d in datasets]
    feature_set = combine_omics_features(dataset_ids)
    pathways = analyze_pathway_convergence(feature_set)
    return CrossOmicsEnrichmentResult(
        program_id=program_id,
        pathways=pathways,
        features=feature_set,
    )

