"""Project export service for selective data export with related entities."""

import hashlib
import io
import json
import zipfile
from datetime import datetime, timezone
from typing import Dict, List, Optional, Set
from uuid import UUID

import pandas as pd
from sqlalchemy.orm import Session

from amprenta_rag.database.models import (
    Compound,
    Dataset,
    Experiment,
    Feature,
    Program,
    Signature,
    SignatureComponent,
)


class ProjectExporter:
    """Service for exporting project data with related entities."""
    
    def __init__(self, db: Session):
        self.db = db
    
    def export_project(
        self,
        program_ids: Optional[List[UUID]] = None,
        experiment_ids: Optional[List[UUID]] = None,
        compound_ids: Optional[List[UUID]] = None,
        include_related: bool = True,
    ) -> bytes:
        """Export project data with related entities to ZIP package.
        
        Args:
            program_ids: List of program UUIDs to export
            experiment_ids: List of experiment UUIDs to export  
            compound_ids: List of compound UUIDs to export
            include_related: Whether to include related entities
            
        Returns:
            ZIP file bytes containing exported data and manifest
        """
        # Collect all entities to export
        export_data = self._collect_export_data(
            program_ids or [],
            experiment_ids or [],
            compound_ids or [],
            include_related
        )
        
        # Create ZIP package
        return self._create_zip_package(export_data)
    
    def _collect_export_data(
        self,
        program_ids: List[UUID],
        experiment_ids: List[UUID],
        compound_ids: List[UUID],
        include_related: bool,
    ) -> Dict[str, List[Dict]]:
        """Collect all data for export."""
        export_data = {
            "programs": [],
            "experiments": [],
            "compounds": [],
            "datasets": [],
            "features": [],
            "signatures": [],
            "signature_components": [],
        }
        
        # Track collected IDs to avoid duplicates
        collected_programs: Set[UUID] = set()
        collected_experiments: Set[UUID] = set()
        collected_compounds: Set[UUID] = set()
        collected_datasets: Set[UUID] = set()
        
        # Export programs
        for program_id in program_ids:
            program = self.db.query(Program).filter(Program.id == program_id).first()
            if program:
                export_data["programs"].append(self._serialize_program(program))
                collected_programs.add(program.id)
                
                if include_related:
                    # Add related experiments
                    for exp in program.experiments:
                        if exp.id not in collected_experiments:
                            export_data["experiments"].append(self._serialize_experiment(exp))
                            collected_experiments.add(exp.id)
                            
                            # Add related datasets
                            for dataset in exp.datasets:
                                if dataset.id not in collected_datasets:
                                    export_data["datasets"].append(self._serialize_dataset(dataset))
                                    collected_datasets.add(dataset.id)
        
        # Export experiments
        for experiment_id in experiment_ids:
            if experiment_id not in collected_experiments:
                experiment = self.db.query(Experiment).filter(Experiment.id == experiment_id).first()
                if experiment:
                    export_data["experiments"].append(self._serialize_experiment(experiment))
                    collected_experiments.add(experiment.id)
                    
                    if include_related:
                        # Add parent program
                        if experiment.program and experiment.program.id not in collected_programs:
                            export_data["programs"].append(self._serialize_program(experiment.program))
                            collected_programs.add(experiment.program.id)
                        
                        # Add related datasets
                        for dataset in experiment.datasets:
                            if dataset.id not in collected_datasets:
                                export_data["datasets"].append(self._serialize_dataset(dataset))
                                collected_datasets.add(dataset.id)
        
        # Export compounds
        for compound_id in compound_ids:
            if compound_id not in collected_compounds:
                compound = self.db.query(Compound).filter(Compound.id == compound_id).first()
                if compound:
                    export_data["compounds"].append(self._serialize_compound(compound))
                    collected_compounds.add(compound.id)
                    
                    if include_related:
                        # Add related features
                        features = self.db.query(Feature).filter(
                            Feature.compound_id == compound.id
                        ).all()
                        for feature in features:
                            export_data["features"].append(self._serialize_feature(feature))
                        
                        # Add related signatures
                        signatures = self.db.query(Signature).filter(
                            Signature.compound_id == compound.id
                        ).all()
                        for signature in signatures:
                            export_data["signatures"].append(self._serialize_signature(signature))
                            
                            # Add signature components
                            components = self.db.query(SignatureComponent).filter(
                                SignatureComponent.signature_id == signature.id
                            ).all()
                            for component in components:
                                export_data["signature_components"].append(
                                    self._serialize_signature_component(component)
                                )
        
        return export_data
    
    def _serialize_program(self, program: Program) -> Dict:
        """Serialize program to dict."""
        return {
            "id": str(program.id),
            "name": program.name,
            "description": program.description,
            "status": program.status,
            "created_at": program.created_at.isoformat() if program.created_at else None,
            "updated_at": program.updated_at.isoformat() if program.updated_at else None,
        }
    
    def _serialize_experiment(self, experiment: Experiment) -> Dict:
        """Serialize experiment to dict."""
        return {
            "id": str(experiment.id),
            "name": experiment.name,
            "description": experiment.description,
            "program_id": str(experiment.program_id) if experiment.program_id else None,
            "experiment_type": experiment.experiment_type,
            "status": experiment.status,
            "created_at": experiment.created_at.isoformat() if experiment.created_at else None,
            "updated_at": experiment.updated_at.isoformat() if experiment.updated_at else None,
        }
    
    def _serialize_compound(self, compound: Compound) -> Dict:
        """Serialize compound to dict."""
        return {
            "id": str(compound.id),
            "smiles": compound.smiles,
            "name": compound.name,
            "molecular_weight": compound.molecular_weight,
            "logp": compound.logp,
            "tpsa": compound.tpsa,
            "hbd": compound.hbd,
            "hba": compound.hba,
            "rotatable_bonds": compound.rotatable_bonds,
            "created_at": compound.created_at.isoformat() if compound.created_at else None,
            "updated_at": compound.updated_at.isoformat() if compound.updated_at else None,
        }
    
    def _serialize_dataset(self, dataset: Dataset) -> Dict:
        """Serialize dataset to dict."""
        return {
            "id": str(dataset.id),
            "name": dataset.name,
            "description": dataset.description,
            "experiment_id": str(dataset.experiment_id) if dataset.experiment_id else None,
            "dataset_type": dataset.dataset_type,
            "file_path": dataset.file_path,
            "file_size": dataset.file_size,
            "row_count": dataset.row_count,
            "column_count": dataset.column_count,
            "created_at": dataset.created_at.isoformat() if dataset.created_at else None,
            "updated_at": dataset.updated_at.isoformat() if dataset.updated_at else None,
        }
    
    def _serialize_feature(self, feature: Feature) -> Dict:
        """Serialize feature to dict."""
        return {
            "id": str(feature.id),
            "name": feature.name,
            "feature_type": feature.feature_type,
            "compound_id": str(feature.compound_id) if feature.compound_id else None,
            "dataset_id": str(feature.dataset_id) if feature.dataset_id else None,
            "value": feature.value,
            "unit": feature.unit,
            "confidence": feature.confidence,
            "created_at": feature.created_at.isoformat() if feature.created_at else None,
        }
    
    def _serialize_signature(self, signature: Signature) -> Dict:
        """Serialize signature to dict."""
        return {
            "id": str(signature.id),
            "name": signature.name,
            "signature_type": signature.signature_type,
            "compound_id": str(signature.compound_id) if signature.compound_id else None,
            "dataset_id": str(signature.dataset_id) if signature.dataset_id else None,
            "metadata": signature.metadata,
            "created_at": signature.created_at.isoformat() if signature.created_at else None,
        }
    
    def _serialize_signature_component(self, component: SignatureComponent) -> Dict:
        """Serialize signature component to dict."""
        return {
            "id": str(component.id),
            "signature_id": str(component.signature_id),
            "feature_name": component.feature_name,
            "value": component.value,
            "rank": component.rank,
            "p_value": component.p_value,
            "fold_change": component.fold_change,
        }
    
    def _create_zip_package(self, export_data: Dict[str, List[Dict]]) -> bytes:
        """Create ZIP package with exported data and manifest."""
        zip_buffer = io.BytesIO()
        
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            # Add data files
            for entity_type, entities in export_data.items():
                if entities:  # Only add non-empty datasets
                    # Convert to DataFrame and CSV
                    df = pd.DataFrame(entities)
                    csv_buffer = io.StringIO()
                    df.to_csv(csv_buffer, index=False)
                    
                    # Add to ZIP
                    zip_file.writestr(
                        f"{entity_type}.csv",
                        csv_buffer.getvalue().encode('utf-8')
                    )
            
            # Create manifest
            manifest = self._create_manifest(export_data, zip_file)
            
            # Add manifest to ZIP
            zip_file.writestr(
                "manifest.json",
                json.dumps(manifest, indent=2).encode('utf-8')
            )
        
        zip_buffer.seek(0)
        return zip_buffer.getvalue()
    
    def _create_manifest(self, export_data: Dict[str, List[Dict]], zip_file: zipfile.ZipFile) -> Dict:
        """Create manifest with export metadata and checksums."""
        manifest = {
            "export_timestamp": datetime.now(timezone.utc).isoformat(),
            "export_type": "project_export",
            "files": {},
            "summary": {},
        }
        
        # Calculate file checksums and sizes
        for info in zip_file.infolist():
            if info.filename != "manifest.json":
                # Read file data to calculate checksum
                file_data = zip_file.read(info.filename)
                checksum = hashlib.sha256(file_data).hexdigest()
                
                manifest["files"][info.filename] = {
                    "size_bytes": info.file_size,
                    "compressed_size_bytes": info.compress_size,
                    "sha256_checksum": checksum,
                }
        
        # Add summary statistics
        for entity_type, entities in export_data.items():
            manifest["summary"][entity_type] = len(entities)
        
        manifest["summary"]["total_records"] = sum(manifest["summary"].values())
        
        return manifest


def export_project(
    program_ids: Optional[List[UUID]] = None,
    experiment_ids: Optional[List[UUID]] = None,
    compound_ids: Optional[List[UUID]] = None,
    db: Optional[Session] = None,
    include_related: bool = True,
) -> bytes:
    """Export project data to ZIP package.
    
    Args:
        program_ids: List of program UUIDs to export
        experiment_ids: List of experiment UUIDs to export
        compound_ids: List of compound UUIDs to export
        db: Database session (required)
        include_related: Whether to include related entities
        
    Returns:
        ZIP file bytes containing exported data and manifest
        
    Raises:
        ValueError: If no database session provided
    """
    if db is None:
        raise ValueError("Database session is required")
    
    exporter = ProjectExporter(db)
    return exporter.export_project(
        program_ids=program_ids,
        experiment_ids=experiment_ids,
        compound_ids=compound_ids,
        include_related=include_related,
    )
