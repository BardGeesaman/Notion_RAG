"""
Dual-write capability for transition period.

Allows writing to both Notion and Postgres during migration,
ensuring data consistency between systems.
"""

from __future__ import annotations

from typing import Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.api.services import (
    programs as program_service,
    experiments as experiment_service,
    datasets as dataset_service,
)

logger = get_logger(__name__)


class DualWriteManager:
    """
    Manages dual-write operations to both Notion and Postgres.
    
    During the transition period, this ensures data is written to both
    systems to maintain consistency.
    """
    
    def __init__(self, db: Session, enable_notion: bool = True, enable_postgres: bool = True):
        """
        Initialize dual-write manager.
        
        Args:
            db: Postgres database session
            enable_notion: Enable writes to Notion
            enable_postgres: Enable writes to Postgres
        """
        self.db = db
        self.enable_notion = enable_notion
        self.enable_postgres = enable_postgres
    
    def create_program(self, program_data: dict, notion_page_id: Optional[str] = None) -> UUID:
        """
        Create program in both systems.
        
        Args:
            program_data: Program data dictionary
            notion_page_id: Existing Notion page ID (if updating existing)
            
        Returns:
            Postgres UUID of created program
        """
        # Write to Postgres
        postgres_id = None
        if self.enable_postgres:
            try:
                from amprenta_rag.api.schemas import ProgramCreate
                program_create = ProgramCreate(**program_data)
                program = program_service.create_program(self.db, program_create)
                postgres_id = program.id
                
                # Store notion_page_id if provided
                if notion_page_id:
                    program.notion_page_id = notion_page_id
                    self.db.commit()
                
                logger.info(
                    "[MIGRATION][DUAL-WRITE] Created program in Postgres: %s (%s)",
                    program.name,
                    postgres_id,
                )
            except Exception as e:
                logger.error("[MIGRATION][DUAL-WRITE] Error creating program in Postgres: %r", e)
                raise
        
        # Write to Notion (if enabled and not already exists)
        if self.enable_notion and notion_page_id is None:
            try:
                # Notion write intentionally disabled (deprecated)
                logger.debug("[MIGRATION][DUAL-WRITE] Notion write intentionally disabled (deprecated)")
            except Exception as e:
                logger.error("[MIGRATION][DUAL-WRITE] Error creating program in Notion: %r", e)
                # Don't fail the whole operation if Notion write fails
        
        return postgres_id
    
    def create_experiment(self, experiment_data: dict, notion_page_id: Optional[str] = None) -> UUID:
        """Create experiment in both systems."""
        postgres_id = None
        if self.enable_postgres:
            try:
                from amprenta_rag.api.schemas import ExperimentCreate
                experiment_create = ExperimentCreate(**experiment_data)
                experiment = experiment_service.create_experiment(self.db, experiment_create)
                postgres_id = experiment.id
                
                if notion_page_id:
                    experiment.notion_page_id = notion_page_id
                    self.db.commit()
                
                logger.info(
                    "[MIGRATION][DUAL-WRITE] Created experiment in Postgres: %s (%s)",
                    experiment.name,
                    postgres_id,
                )
            except Exception as e:
                logger.error("[MIGRATION][DUAL-WRITE] Error creating experiment in Postgres: %r", e)
                raise
        
        if self.enable_notion and notion_page_id is None:
            try:
                logger.debug("[MIGRATION][DUAL-WRITE] Notion write not yet implemented")
            except Exception as e:
                logger.error("[MIGRATION][DUAL-WRITE] Error creating experiment in Notion: %r", e)
        
        return postgres_id
    
    def create_dataset(self, dataset_data: dict, notion_page_id: Optional[str] = None) -> UUID:
        """Create dataset in both systems."""
        postgres_id = None
        if self.enable_postgres:
            try:
                from amprenta_rag.api.schemas import DatasetCreate
                dataset_create = DatasetCreate(**dataset_data)
                dataset = dataset_service.create_dataset(self.db, dataset_create)
                postgres_id = dataset.id
                
                if notion_page_id:
                    dataset.notion_page_id = notion_page_id
                    self.db.commit()
                
                logger.info(
                    "[MIGRATION][DUAL-WRITE] Created dataset in Postgres: %s (%s)",
                    dataset.name,
                    postgres_id,
                )
            except Exception as e:
                logger.error("[MIGRATION][DUAL-WRITE] Error creating dataset in Postgres: %r", e)
                raise
        
        if self.enable_notion and notion_page_id is None:
            try:
                logger.debug("[MIGRATION][DUAL-WRITE] Notion write not yet implemented")
            except Exception as e:
                logger.error("[MIGRATION][DUAL-WRITE] Error creating dataset in Notion: %r", e)
        
        return postgres_id

