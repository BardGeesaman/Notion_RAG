"""Data lineage visualization utilities."""
from __future__ import annotations

from typing import Dict, List, Optional, Set
from uuid import UUID

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _add_node(
    nodes: List[Dict],
    node_id: str,
    node_type: str,
    name: str,
    metadata: Optional[Dict] = None,
) -> None:
    """
    Add a node to the nodes list if it doesn't already exist.

    Args:
        nodes: List of node dictionaries
        node_id: Unique identifier for the node
        node_type: Type of entity (experiment, dataset, signature, etc.)
        name: Display name
        metadata: Optional metadata dictionary
    """
    # Check if node already exists
    if any(node["id"] == node_id for node in nodes):
        return

    nodes.append({
        "id": node_id,
        "type": node_type,
        "name": name,
        "metadata": metadata or {},
    })


def _add_edge(
    edges: List[Dict],
    source: str,
    target: str,
    relationship: str,
) -> None:
    """
    Add an edge to the edges list if it doesn't already exist.

    Args:
        edges: List of edge dictionaries
        source: Source node ID
        target: Target node ID
        relationship: Type of relationship
    """
    # Check if edge already exists
    if any(
        edge["source"] == source and edge["target"] == target and edge["relationship"] == relationship
        for edge in edges
    ):
        return

    edges.append({
        "source": source,
        "target": target,
        "relationship": relationship,
    })


def get_experiment_lineage(experiment_id: UUID, db) -> Dict[str, List[Dict]]:
    """
    Get lineage graph for an experiment.

    Args:
        experiment_id: UUID of the experiment
        db: Database session

    Returns:
        Dictionary with "nodes" and "edges" lists
    """
    from amprenta_rag.database.models import Experiment, Project

    nodes: List[Dict] = []
    edges: List[Dict] = []

    # Load experiment
    experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()
    if not experiment:
        return {"nodes": [], "edges": []}

    exp_id_str = str(experiment.id)

    # Add experiment as root node
    _add_node(
        nodes,
        exp_id_str,
        "experiment",
        experiment.name or "Unknown Experiment",
        {
            "design_type": experiment.design_type,
            "created_at": experiment.created_at.isoformat() if experiment.created_at else None,
        },
    )

    # Add connected datasets
    if hasattr(experiment, "datasets") and experiment.datasets:
        for dataset in experiment.datasets:
            dataset_id_str = str(dataset.id)
            _add_node(
                nodes,
                dataset_id_str,
                "dataset",
                dataset.name or "Unknown Dataset",
                {
                    "omics_type": dataset.omics_type,
                    "created_at": dataset.created_at.isoformat() if dataset.created_at else None,
                },
            )
            _add_edge(edges, exp_id_str, dataset_id_str, "has_dataset")

            # Add signatures connected to this dataset
            if hasattr(dataset, "signatures") and dataset.signatures:
                for signature in dataset.signatures:
                    sig_id_str = str(signature.id)
                    _add_node(
                        nodes,
                        sig_id_str,
                        "signature",
                        signature.name or "Unknown Signature",
                        {
                            "created_at": signature.created_at.isoformat() if signature.created_at else None,
                        },
                    )
                    _add_edge(edges, dataset_id_str, sig_id_str, "matches_signature")

    # Add linked protocols
    if hasattr(experiment, "protocol_links") and experiment.protocol_links:
        for protocol_link in experiment.protocol_links:
            if hasattr(protocol_link, "protocol") and protocol_link.protocol:
                protocol = protocol_link.protocol
                protocol_id_str = str(protocol.id)
                _add_node(
                    nodes,
                    protocol_id_str,
                    "protocol",
                    protocol.name or "Unknown Protocol",
                    {
                        "category": protocol.category,
                        "version": protocol.version,
                    },
                )
                _add_edge(edges, exp_id_str, protocol_id_str, "uses_protocol")

    # Add project if linked
    if experiment.project_id:
        project = db.query(Project).filter(Project.id == experiment.project_id).first()
        if project:
            project_id_str = str(project.id)
            _add_node(
                nodes,
                project_id_str,
                "project",
                project.name or "Unknown Project",
                {
                    "created_at": project.created_at.isoformat() if project.created_at else None,
                },
            )
            _add_edge(edges, exp_id_str, project_id_str, "belongs_to_project")

    return {"nodes": nodes, "edges": edges}


def get_entity_lineage(entity_type: str, entity_id: UUID, db) -> Dict[str, List[Dict]]:
    """
    Get lineage graph for any entity type.

    Args:
        entity_type: Type of entity (experiment, dataset, signature, etc.)
        entity_id: UUID of the entity
        db: Database session

    Returns:
        Dictionary with "nodes" and "edges" lists
    """
    if entity_type.lower() == "experiment":
        return get_experiment_lineage(entity_id, db)
    elif entity_type.lower() == "dataset":
        return get_dataset_lineage(entity_id, db)
    elif entity_type.lower() == "signature":
        return get_signature_lineage(entity_id, db)
    else:
        logger.warning(f"Unsupported entity type for lineage: {entity_type}")
        return {"nodes": [], "edges": []}


def get_dataset_lineage(dataset_id: UUID, db) -> Dict[str, List[Dict]]:
    """
    Get lineage graph for a dataset.

    Args:
        dataset_id: UUID of the dataset
        db: Database session

    Returns:
        Dictionary with "nodes" and "edges" lists
    """
    from amprenta_rag.database.models import Dataset

    nodes: List[Dict] = []
    edges: List[Dict] = []

    # Load dataset
    dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
    if not dataset:
        return {"nodes": [], "edges": []}

    dataset_id_str = str(dataset.id)

    # Add dataset as root node
    _add_node(
        nodes,
        dataset_id_str,
        "dataset",
        dataset.name or "Unknown Dataset",
        {
            "omics_type": dataset.omics_type,
            "created_at": dataset.created_at.isoformat() if dataset.created_at else None,
        },
    )

    # Add connected experiments
    if hasattr(dataset, "experiments") and dataset.experiments:
        for experiment in dataset.experiments:
            exp_id_str = str(experiment.id)
            _add_node(
                nodes,
                exp_id_str,
                "experiment",
                experiment.name or "Unknown Experiment",
                {
                    "design_type": experiment.design_type,
                },
            )
            _add_edge(edges, exp_id_str, dataset_id_str, "has_dataset")

    # Add connected signatures
    if hasattr(dataset, "signatures") and dataset.signatures:
        for signature in dataset.signatures:
            sig_id_str = str(signature.id)
            _add_node(
                nodes,
                sig_id_str,
                "signature",
                signature.name or "Unknown Signature",
                {},
            )
            _add_edge(edges, dataset_id_str, sig_id_str, "matches_signature")

    return {"nodes": nodes, "edges": edges}


def get_signature_lineage(signature_id: UUID, db) -> Dict[str, List[Dict]]:
    """
    Get lineage graph for a signature.

    Args:
        signature_id: UUID of the signature
        db: Database session

    Returns:
        Dictionary with "nodes" and "edges" lists
    """
    from amprenta_rag.database.models import Signature

    nodes: List[Dict] = []
    edges: List[Dict] = []

    # Load signature
    signature = db.query(Signature).filter(Signature.id == signature_id).first()
    if not signature:
        return {"nodes": [], "edges": []}

    sig_id_str = str(signature.id)

    # Add signature as root node
    _add_node(
        nodes,
        sig_id_str,
        "signature",
        signature.name or "Unknown Signature",
        {
            "created_at": signature.created_at.isoformat() if signature.created_at else None,
        },
    )

    # Add connected datasets
    if hasattr(signature, "datasets") and signature.datasets:
        for dataset in signature.datasets:
            dataset_id_str = str(dataset.id)
            _add_node(
                nodes,
                dataset_id_str,
                "dataset",
                dataset.name or "Unknown Dataset",
                {
                    "omics_type": dataset.omics_type,
                },
            )
            _add_edge(edges, dataset_id_str, sig_id_str, "matches_signature")

            # Add experiments connected to datasets
            if hasattr(dataset, "experiments") and dataset.experiments:
                for experiment in dataset.experiments:
                    exp_id_str = str(experiment.id)
                    _add_node(
                        nodes,
                        exp_id_str,
                        "experiment",
                        experiment.name or "Unknown Experiment",
                        {},
                    )
                    _add_edge(edges, exp_id_str, dataset_id_str, "has_dataset")

    return {"nodes": nodes, "edges": edges}


def get_full_lineage_graph(db, limit: int = 100) -> Dict[str, List[Dict]]:
    """
    Get a combined lineage graph of recent entities.

    Args:
        db: Database session
        limit: Maximum number of experiments to include

    Returns:
        Dictionary with "nodes" and "edges" lists
    """
    from amprenta_rag.database.models import Experiment

    nodes: List[Dict] = []
    edges: List[Dict] = []
    seen_nodes: Set[str] = set()
    seen_edges: Set[str] = set()

    # Get recent experiments
    experiments = db.query(Experiment).order_by(Experiment.created_at.desc()).limit(limit).all()

    for experiment in experiments:
        # Get lineage for this experiment
        exp_lineage = get_experiment_lineage(experiment.id, db)

        # Merge nodes (avoid duplicates)
        for node in exp_lineage["nodes"]:
            node_id = node["id"]
            if node_id not in seen_nodes:
                nodes.append(node)
                seen_nodes.add(node_id)

        # Merge edges (avoid duplicates)
        for edge in exp_lineage["edges"]:
            edge_key = f"{edge['source']}-{edge['target']}-{edge['relationship']}"
            if edge_key not in seen_edges:
                edges.append(edge)
                seen_edges.add(edge_key)

    return {"nodes": nodes, "edges": edges}

