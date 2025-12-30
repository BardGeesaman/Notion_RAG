"""HTS plate and well integration for imaging data."""

from __future__ import annotations

import logging
import re
from typing import List, Optional, Dict, Any
from uuid import UUID

from sqlalchemy.orm import Session
from sqlalchemy import and_, func

from amprenta_rag.models.chemistry import HTSCampaign, HTSPlate, HTSWell, Compound
from amprenta_rag.imaging.models import MicroscopyImage, CellSegmentation, CellFeature
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


class HTSImagingIntegration:
    """Integration between HTS plate/well data and imaging workflow."""
    
    def __init__(self):
        """Initialize HTS imaging integration."""
        pass
    
    def link_images_to_plate(
        self,
        plate_barcode: str,
        images: List[Dict[str, Any]],
        campaign_id: UUID,
        db: Session,
        auto_create_wells: bool = True
    ) -> Dict[str, Any]:
        """
        Link microscopy images to HTS plate wells.
        
        Args:
            plate_barcode: Unique plate identifier
            images: List of image metadata with well positions
            campaign_id: HTS campaign UUID
            db: Database session
            auto_create_wells: Whether to auto-create missing wells
        
        Returns:
            Dictionary with linking results and statistics
        """
        try:
            logger.info(f"Linking {len(images)} images to plate {plate_barcode}")
            
            # Get or create plate
            plate = db.query(HTSPlate).filter(HTSPlate.barcode == plate_barcode).first()
            if not plate:
                # Create new plate
                plate = HTSPlate(
                    campaign_id=campaign_id,
                    barcode=plate_barcode,
                    plate_format="384",  # Default format
                    layout_metadata={}
                )
                db.add(plate)
                db.commit()
                db.refresh(plate)
                logger.info(f"Created new plate {plate_barcode}")
            
            # Track results
            results = {
                "plate_id": plate.id,
                "plate_barcode": plate_barcode,
                "images_processed": 0,
                "images_linked": 0,
                "wells_created": 0,
                "wells_found": 0,
                "errors": []
            }
            
            # Process each image
            for image_data in images:
                try:
                    well_position = image_data.get("well_position", "").upper()
                    image_path = image_data.get("image_path", "")
                    channel = image_data.get("channel", "")
                    
                    if not well_position or not image_path or not channel:
                        error_msg = f"Missing required fields in image data: {image_data}"
                        logger.warning(error_msg)
                        results["errors"].append(error_msg)
                        continue
                    
                    # Validate well position format (e.g., A01, P24)
                    if not self._is_valid_well_position(well_position):
                        error_msg = f"Invalid well position format: {well_position}"
                        logger.warning(error_msg)
                        results["errors"].append(error_msg)
                        continue
                    
                    # Get or create well
                    well = db.query(HTSWell).filter(
                        and_(
                            HTSWell.plate_id == plate.id,
                            HTSWell.position == well_position
                        )
                    ).first()
                    
                    if not well and auto_create_wells:
                        well = HTSWell(
                            plate_id=plate.id,
                            position=well_position,
                            compound_id=image_data.get("compound_id"),
                            concentration=image_data.get("concentration"),
                            treatment_metadata=image_data.get("treatment_metadata", {})
                        )
                        db.add(well)
                        db.commit()
                        db.refresh(well)
                        results["wells_created"] += 1
                        logger.debug(f"Created well {well_position} in plate {plate_barcode}")
                    elif well:
                        results["wells_found"] += 1
                    else:
                        error_msg = f"Well {well_position} not found and auto_create_wells=False"
                        logger.warning(error_msg)
                        results["errors"].append(error_msg)
                        continue
                    
                    # Check if image already exists
                    existing_image = db.query(MicroscopyImage).filter(
                        and_(
                            MicroscopyImage.well_id == well.id,
                            MicroscopyImage.channel == channel,
                            MicroscopyImage.z_slice == image_data.get("z_slice", 0),
                            MicroscopyImage.timepoint == image_data.get("timepoint", 0)
                        )
                    ).first()
                    
                    if existing_image:
                        logger.debug(f"Image already exists for well {well_position}, channel {channel}")
                        results["images_processed"] += 1
                        continue
                    
                    # Create microscopy image record
                    image = MicroscopyImage(
                        well_id=well.id,
                        channel=channel,
                        z_slice=image_data.get("z_slice", 0),
                        timepoint=image_data.get("timepoint", 0),
                        width=image_data.get("width", 0),
                        height=image_data.get("height", 0),
                        bit_depth=image_data.get("bit_depth", 16),
                        pixel_size_um=image_data.get("pixel_size_um"),
                        image_path=image_path,
                        thumbnail_path=image_data.get("thumbnail_path"),
                        image_metadata=image_data.get("metadata", {}),
                        focus_score=image_data.get("focus_score"),
                        signal_noise_ratio=image_data.get("signal_noise_ratio"),
                        acquired_at=image_data.get("acquired_at")
                    )
                    
                    db.add(image)
                    results["images_linked"] += 1
                    results["images_processed"] += 1
                    
                    logger.debug(f"Linked image {image_path} to well {well_position}")
                    
                except Exception as e:
                    error_msg = f"Failed to process image {image_data}: {str(e)}"
                    logger.error(error_msg)
                    results["errors"].append(error_msg)
                    continue
            
            # Commit all changes
            db.commit()
            
            logger.info(
                f"Plate linking complete: {results['images_linked']}/{results['images_processed']} "
                f"images linked, {results['wells_created']} wells created"
            )
            
            return results
            
        except Exception as e:
            logger.error(f"Failed to link images to plate {plate_barcode}: {str(e)}")
            db.rollback()
            raise
    
    def get_plate_images(
        self,
        plate_id: UUID,
        db: Session,
        channel: Optional[str] = None,
        z_slice: Optional[int] = None,
        timepoint: Optional[int] = None
    ) -> List[MicroscopyImage]:
        """
        Get all microscopy images for a plate.
        
        Args:
            plate_id: Plate UUID
            db: Database session
            channel: Optional channel filter
            z_slice: Optional z-slice filter
            timepoint: Optional timepoint filter
        
        Returns:
            List of microscopy images
        """
        try:
            # Build query
            query = db.query(MicroscopyImage).join(HTSWell).filter(
                HTSWell.plate_id == plate_id
            )
            
            # Apply filters
            if channel:
                query = query.filter(MicroscopyImage.channel == channel)
            if z_slice is not None:
                query = query.filter(MicroscopyImage.z_slice == z_slice)
            if timepoint is not None:
                query = query.filter(MicroscopyImage.timepoint == timepoint)
            
            # Order by well position and channel
            images = query.order_by(HTSWell.position, MicroscopyImage.channel).all()
            
            logger.info(f"Retrieved {len(images)} images for plate {plate_id}")
            return images
            
        except Exception as e:
            logger.error(f"Failed to get plate images: {str(e)}")
            raise
    
    def get_well_images(
        self,
        well_id: UUID,
        db: Session,
        channel: Optional[str] = None
    ) -> List[MicroscopyImage]:
        """
        Get all microscopy images for a specific well.
        
        Args:
            well_id: Well UUID
            db: Database session
            channel: Optional channel filter
        
        Returns:
            List of microscopy images
        """
        try:
            query = db.query(MicroscopyImage).filter(MicroscopyImage.well_id == well_id)
            
            if channel:
                query = query.filter(MicroscopyImage.channel == channel)
            
            images = query.order_by(
                MicroscopyImage.channel,
                MicroscopyImage.z_slice,
                MicroscopyImage.timepoint
            ).all()
            
            logger.debug(f"Retrieved {len(images)} images for well {well_id}")
            return images
            
        except Exception as e:
            logger.error(f"Failed to get well images: {str(e)}")
            raise
    
    def get_plate_well_positions(self, plate_id: UUID, db: Session) -> List[str]:
        """
        Get all well positions for a plate.
        
        Args:
            plate_id: Plate UUID
            db: Database session
        
        Returns:
            List of well position strings (e.g., ["A01", "A02", ...])
        """
        try:
            positions = db.query(HTSWell.position).filter(
                HTSWell.plate_id == plate_id
            ).order_by(HTSWell.position).all()
            
            return [pos[0] for pos in positions]
            
        except Exception as e:
            logger.error(f"Failed to get plate well positions: {str(e)}")
            raise
    
    def get_plate_channels(self, plate_id: UUID, db: Session) -> List[str]:
        """
        Get all available channels for a plate.
        
        Args:
            plate_id: Plate UUID
            db: Database session
        
        Returns:
            List of channel names
        """
        try:
            channels = db.query(MicroscopyImage.channel).join(HTSWell).filter(
                HTSWell.plate_id == plate_id
            ).distinct().order_by(MicroscopyImage.channel).all()
            
            return [ch[0] for ch in channels]
            
        except Exception as e:
            logger.error(f"Failed to get plate channels: {str(e)}")
            raise
    
    def get_plate_statistics(self, plate_id: UUID, db: Session) -> Dict[str, Any]:
        """
        Get imaging statistics for a plate.
        
        Args:
            plate_id: Plate UUID
            db: Database session
        
        Returns:
            Dictionary with plate imaging statistics
        """
        try:
            # Get plate info
            plate = db.query(HTSPlate).filter(HTSPlate.id == plate_id).first()
            if not plate:
                raise ValueError(f"Plate {plate_id} not found")
            
            # Count wells and images
            well_count = db.query(HTSWell).filter(HTSWell.plate_id == plate_id).count()
            
            image_count = db.query(MicroscopyImage).join(HTSWell).filter(
                HTSWell.plate_id == plate_id
            ).count()
            
            # Count segmentations
            segmentation_count = db.query(CellSegmentation).join(MicroscopyImage).join(HTSWell).filter(
                HTSWell.plate_id == plate_id
            ).count()
            
            # Get channels
            channels = self.get_plate_channels(plate_id, db)
            
            # Calculate imaging coverage
            wells_with_images = db.query(HTSWell.id).join(MicroscopyImage).filter(
                HTSWell.plate_id == plate_id
            ).distinct().count()
            
            imaging_coverage = wells_with_images / well_count if well_count > 0 else 0.0
            
            return {
                "plate_id": plate_id,
                "plate_barcode": plate.barcode,
                "plate_format": plate.plate_format,
                "well_count": well_count,
                "image_count": image_count,
                "segmentation_count": segmentation_count,
                "channels": channels,
                "wells_with_images": wells_with_images,
                "imaging_coverage": imaging_coverage,
                "images_per_well": image_count / well_count if well_count > 0 else 0.0
            }
            
        except Exception as e:
            logger.error(f"Failed to get plate statistics: {str(e)}")
            raise
    
    def _is_valid_well_position(self, position: str) -> bool:
        """
        Validate well position format.
        
        Args:
            position: Well position string
        
        Returns:
            True if valid format (e.g., A01, P24)
        """
        # Pattern: Letter(s) followed by digits (e.g., A01, AA01, P24)
        pattern = r'^[A-Z]+[0-9]+$'
        return bool(re.match(pattern, position.upper()))
    
    def parse_well_position(self, position: str) -> Dict[str, int]:
        """
        Parse well position into row and column indices.
        
        Args:
            position: Well position (e.g., "A01", "P24")
        
        Returns:
            Dictionary with row and column indices (0-based)
        """
        position = position.upper()
        
        # Extract letter part (row) and number part (column)
        match = re.match(r'^([A-Z]+)([0-9]+)$', position)
        if not match:
            raise ValueError(f"Invalid well position format: {position}")
        
        row_letters, col_digits = match.groups()
        
        # Convert letters to row index (A=0, B=1, ..., AA=26, etc.)
        row_index = 0
        for i, letter in enumerate(reversed(row_letters)):
            row_index += (ord(letter) - ord('A') + 1) * (26 ** i)
        row_index -= 1  # Convert to 0-based
        
        # Convert digits to column index (1-based to 0-based)
        col_index = int(col_digits) - 1
        
        return {
            "row": row_index,
            "column": col_index,
            "row_letter": row_letters,
            "column_number": int(col_digits)
        }
