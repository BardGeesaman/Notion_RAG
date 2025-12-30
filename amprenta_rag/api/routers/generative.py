"""Generative chemistry API endpoints."""

from __future__ import annotations

from typing import List

from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.security import HTTPBearer

from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.api.schemas import (
    SampleRequest,
    SampleResponse,
    InterpolateRequest,
    InterpolateResponse,
    OptimizeRequest,
    OptimizeResponse,
    ScaffoldHopRequest,
    ScaffoldHopResponse,
    GenerativeModelsResponse,
    GenerativeModelInfo,
    MoleculeSchema,
    ErrorResponse,
)
from amprenta_rag.ml.generative.service import GenerativeChemistryService, GenerativeModelNotFoundError
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)
security = HTTPBearer()

router = APIRouter()

# Global service instance (would be initialized with model path in production)
_service: GenerativeChemistryService = None


def get_generative_service() -> GenerativeChemistryService:
    """Get or create generative chemistry service instance."""
    global _service
    if _service is None:
        # In production, this would load from configuration
        _service = GenerativeChemistryService(
            model_path=None,  # No model path - will return appropriate errors
        )
    return _service


@router.post(
    "/generative/sample",
    response_model=SampleResponse,
    status_code=status.HTTP_200_OK,
    summary="Generate random molecules",
    description="Generate random molecules by sampling from the VAE latent space."
)
async def sample_molecules(
    request: SampleRequest,
    current_user=Depends(get_current_user),
    service: GenerativeChemistryService = Depends(get_generative_service),
) -> SampleResponse:
    """Generate random molecules from latent space."""
    try:
        logger.info(f"User {current_user} requesting {request.n_samples} random molecules")
        
        # Generate molecules
        molecules_data = service.sample(
            n_samples=request.n_samples,
            temperature=request.temperature,
            max_length=request.max_length,
        )
        
        # Convert to response format
        molecules = [MoleculeSchema(**mol) for mol in molecules_data]
        model_info = service.get_model_info()
        
        return SampleResponse(
            molecules=molecules,
            count=len(molecules),
            model_info=model_info,
        )
    
    except GenerativeModelNotFoundError as e:
        logger.error(f"Model not found: {e}")
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Generative model not available"
        )
    except Exception as e:
        logger.error(f"Failed to generate molecules: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to generate molecules"
        )


@router.post(
    "/generative/interpolate",
    response_model=InterpolateResponse,
    status_code=status.HTTP_200_OK,
    summary="Interpolate between molecules",
    description="Generate molecules by interpolating between two input molecules in latent space."
)
async def interpolate_molecules(
    request: InterpolateRequest,
    current_user=Depends(get_current_user),
    service: GenerativeChemistryService = Depends(get_generative_service),
) -> InterpolateResponse:
    """Interpolate between two molecules in latent space."""
    try:
        logger.info(f"User {current_user} interpolating between {request.smiles_start} and {request.smiles_end}")
        
        # Generate interpolated molecules
        molecules_data = service.interpolate(
            smiles_start=request.smiles_start,
            smiles_end=request.smiles_end,
            steps=request.steps,
            interpolation_type=request.interpolation_type,
        )
        
        # Convert to response format
        molecules = [MoleculeSchema(**mol) for mol in molecules_data]
        
        return InterpolateResponse(
            molecules=molecules,
            start_smiles=request.smiles_start,
            end_smiles=request.smiles_end,
            steps=request.steps,
        )
    
    except GenerativeModelNotFoundError as e:
        logger.error(f"Model not found: {e}")
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Generative model not available"
        )
    except ValueError as e:
        logger.warning(f"Invalid input: {e}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )
    except Exception as e:
        logger.error(f"Failed to interpolate molecules: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to interpolate molecules"
        )


@router.post(
    "/generative/optimize",
    response_model=OptimizeResponse,
    status_code=status.HTTP_200_OK,
    summary="Optimize molecular properties",
    description="Optimize molecules for desired properties using property-guided generation."
)
async def optimize_molecules(
    request: OptimizeRequest,
    current_user=Depends(get_current_user),
    service: GenerativeChemistryService = Depends(get_generative_service),
) -> OptimizeResponse:
    """Optimize molecules for desired properties."""
    try:
        logger.info(f"User {current_user} optimizing from seed {request.seed_smiles}")
        
        # Convert constraints to dict format
        constraints_data = []
        for constraint in request.constraints:
            constraints_data.append({
                "name": constraint.name,
                "min_value": constraint.min_value,
                "max_value": constraint.max_value,
                "target_value": constraint.target_value,
                "weight": constraint.weight,
            })
        
        # Run optimization
        result = service.optimize(
            seed_smiles=request.seed_smiles,
            constraints=constraints_data,
            n_iterations=request.n_iterations,
            n_samples_per_iter=request.n_samples_per_iter,
            learning_rate=request.learning_rate,
            temperature=request.temperature,
        )
        
        # Convert to response format
        optimized_molecules = [MoleculeSchema(**mol) for mol in result["optimized"]]
        
        return OptimizeResponse(
            optimized=optimized_molecules,
            seed_smiles=request.seed_smiles,
            seed_properties=result["seed_properties"],
            best_score=result["best_score"],
            iterations_completed=result["iterations_completed"],
        )
    
    except GenerativeModelNotFoundError as e:
        logger.error(f"Model not found: {e}")
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Generative model not available"
        )
    except ValueError as e:
        logger.warning(f"Invalid input: {e}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )
    except Exception as e:
        logger.error(f"Failed to optimize molecules: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to optimize molecules"
        )


@router.post(
    "/generative/scaffold-hop",
    response_model=ScaffoldHopResponse,
    status_code=status.HTTP_200_OK,
    summary="Generate scaffold analogs",
    description="Generate molecular analogs by scaffold hopping or preservation."
)
async def scaffold_hop(
    request: ScaffoldHopRequest,
    current_user=Depends(get_current_user),
    service: GenerativeChemistryService = Depends(get_generative_service),
) -> ScaffoldHopResponse:
    """Generate scaffold analogs or novel scaffolds."""
    try:
        logger.info(f"User {current_user} scaffold hopping for {request.smiles}")
        
        # Generate scaffold analogs
        result = service.scaffold_hop(
            smiles=request.smiles,
            n_analogs=request.n_analogs,
            preserve_scaffold=request.preserve_scaffold,
            similarity_threshold=request.similarity_threshold,
        )
        
        # Convert to response format
        analogs = [MoleculeSchema(**mol) for mol in result["analogs"]]
        
        return ScaffoldHopResponse(
            scaffold=result["scaffold"],
            analogs=analogs,
            input_smiles=request.smiles,
            n_generated=result["n_generated"],
        )
    
    except GenerativeModelNotFoundError as e:
        logger.error(f"Model not found: {e}")
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Generative model not available"
        )
    except ValueError as e:
        logger.warning(f"Invalid input: {e}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )
    except Exception as e:
        logger.error(f"Failed to generate scaffold analogs: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to generate scaffold analogs"
        )


@router.get(
    "/generative/models",
    response_model=GenerativeModelsResponse,
    status_code=status.HTTP_200_OK,
    summary="List available models",
    description="List all available generative chemistry models."
)
async def list_models(
    current_user=Depends(get_current_user),
    service: GenerativeChemistryService = Depends(get_generative_service),
) -> GenerativeModelsResponse:
    """List available generative chemistry models."""
    try:
        logger.info(f"User {current_user} requesting available models")
        
        # Get model info
        model_info = service.get_model_info()
        
        # Convert to model info schema
        model = GenerativeModelInfo(
            name=model_info.get("name", "default"),
            version=model_info.get("version", "1.0.0"),
            latent_dim=model_info.get("latent_dim", 0),
            vocab_size=model_info.get("vocab_size", 0),
            status=model_info.get("status", "unknown"),
            description="Molecular VAE for de novo drug design",
        )
        
        return GenerativeModelsResponse(
            models=[model],
            count=1,
            default_model=model.name if model.status == "loaded" else None,
        )
    
    except Exception as e:
        logger.error(f"Failed to list models: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to list models"
        )
