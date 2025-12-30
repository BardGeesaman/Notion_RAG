"""Pydantic schemas for generative chemistry API."""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from pydantic import BaseModel, Field, field_validator


class PropertyConstraintSchema(BaseModel):
    """Schema for property constraint in optimization."""
    
    name: str = Field(..., description="Property name (e.g., 'logP', 'herg', 'mw')")
    min_value: Optional[float] = Field(None, description="Minimum allowed value")
    max_value: Optional[float] = Field(None, description="Maximum allowed value")
    target_value: Optional[float] = Field(None, description="Target value for optimization")
    weight: float = Field(1.0, ge=0.0, description="Importance weight")
    
    @field_validator('name')
    @classmethod
    def validate_name(cls, v):
        if not v or not v.strip():
            raise ValueError("Property name cannot be empty")
        return v.strip().lower()


class MoleculeSchema(BaseModel):
    """Schema for generated molecule with properties."""
    
    smiles: str = Field(..., description="Generated SMILES string")
    properties: dict[str, float] = Field(default_factory=dict, description="Predicted properties")
    score: Optional[float] = Field(None, description="Overall optimization score")
    step: Optional[int] = Field(None, description="Generation step (for interpolation)")
    iteration: Optional[int] = Field(None, description="Optimization iteration")


class SampleRequest(BaseModel):
    """Request schema for random sampling."""
    
    n_samples: int = Field(10, ge=1, le=100, description="Number of molecules to generate")
    temperature: float = Field(1.0, ge=0.1, le=2.0, description="Sampling temperature")
    max_length: int = Field(100, ge=20, le=200, description="Maximum SMILES length")


class SampleResponse(BaseModel):
    """Response schema for random sampling."""
    
    molecules: List[MoleculeSchema] = Field(..., description="Generated molecules")
    count: int = Field(..., description="Number of molecules generated")
    model_info: dict = Field(..., description="Model information used")


class InterpolateRequest(BaseModel):
    """Request schema for molecular interpolation."""
    
    smiles_start: str = Field(..., description="Starting molecule SMILES")
    smiles_end: str = Field(..., description="Ending molecule SMILES")
    steps: int = Field(10, ge=2, le=50, description="Number of interpolation steps")
    interpolation_type: str = Field("linear", description="Interpolation method")
    
    @field_validator('smiles_start', 'smiles_end')
    @classmethod
    def validate_smiles(cls, v):
        if not v or not v.strip():
            raise ValueError("SMILES cannot be empty")
        return v.strip()


class InterpolateResponse(BaseModel):
    """Response schema for molecular interpolation."""
    
    molecules: List[MoleculeSchema] = Field(..., description="Interpolated molecules")
    start_smiles: str = Field(..., description="Starting molecule")
    end_smiles: str = Field(..., description="Ending molecule")
    steps: int = Field(..., description="Number of steps")


class OptimizeRequest(BaseModel):
    """Request schema for property-guided optimization."""
    
    seed_smiles: str = Field(..., description="Starting molecule SMILES")
    constraints: List[PropertyConstraintSchema] = Field(..., description="Property constraints")
    n_iterations: int = Field(100, ge=1, le=1000, description="Number of optimization iterations")
    n_samples_per_iter: int = Field(10, ge=1, le=50, description="Samples per iteration")
    learning_rate: float = Field(0.1, ge=0.01, le=1.0, description="Optimization learning rate")
    temperature: float = Field(1.0, ge=0.1, le=2.0, description="Initial sampling temperature")
    
    @field_validator('seed_smiles')
    @classmethod
    def validate_seed_smiles(cls, v):
        if not v or not v.strip():
            raise ValueError("Seed SMILES cannot be empty")
        return v.strip()
    
    @field_validator('constraints')
    @classmethod
    def validate_constraints(cls, v):
        if not v:
            raise ValueError("At least one constraint must be specified")
        return v


class OptimizeResponse(BaseModel):
    """Response schema for property-guided optimization."""
    
    optimized: List[MoleculeSchema] = Field(..., description="Optimized molecules")
    seed_smiles: str = Field(..., description="Starting molecule")
    seed_properties: dict[str, float] = Field(..., description="Properties of seed molecule")
    best_score: float = Field(..., description="Best optimization score achieved")
    iterations_completed: int = Field(..., description="Number of iterations completed")


class ScaffoldHopRequest(BaseModel):
    """Request schema for scaffold hopping."""
    
    smiles: str = Field(..., description="Input molecule SMILES")
    n_analogs: int = Field(20, ge=1, le=50, description="Number of analogs to generate")
    preserve_scaffold: bool = Field(True, description="Whether to preserve the scaffold")
    similarity_threshold: float = Field(0.7, ge=0.1, le=0.9, description="Similarity threshold")
    
    @field_validator('smiles')
    @classmethod
    def validate_smiles(cls, v):
        if not v or not v.strip():
            raise ValueError("SMILES cannot be empty")
        return v.strip()


class ScaffoldHopResponse(BaseModel):
    """Response schema for scaffold hopping."""
    
    scaffold: Optional[str] = Field(None, description="Extracted scaffold SMILES")
    analogs: List[MoleculeSchema] = Field(..., description="Generated analogs")
    input_smiles: str = Field(..., description="Original input molecule")
    n_generated: int = Field(..., description="Number of analogs generated")


class GenerativeModelInfo(BaseModel):
    """Schema for generative model information."""
    
    name: str = Field(..., description="Model name")
    version: str = Field(..., description="Model version")
    latent_dim: int = Field(..., description="Latent space dimensionality")
    vocab_size: int = Field(..., description="Vocabulary size")
    status: str = Field(..., description="Model status (loaded/not_found/error)")
    created_at: Optional[str] = Field(None, description="Model creation timestamp")
    description: Optional[str] = Field(None, description="Model description")


class GenerativeModelsResponse(BaseModel):
    """Response schema for listing available models."""
    
    models: List[GenerativeModelInfo] = Field(..., description="Available generative models")
    count: int = Field(..., description="Number of models")
    default_model: Optional[str] = Field(None, description="Default model name")


class ErrorResponse(BaseModel):
    """Schema for error responses."""
    
    error: str = Field(..., description="Error message")
    detail: Optional[str] = Field(None, description="Detailed error information")
    error_type: str = Field(..., description="Type of error")
    request_id: Optional[str] = Field(None, description="Request identifier for tracking")
