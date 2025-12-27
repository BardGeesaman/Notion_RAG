"""
FastAPI main application.

This module sets up the FastAPI application with all routes, middleware,
and configuration.
"""

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from amprenta_rag.api.routers import (
    datasets,
    experiments,
    features,
    programs,
    signatures,
    compounds,
    screening,
    sar,
    reports,
    explainability,
    analysis,
    batch,
    quality,
    protocols,
    hts,
    pathways,
    moa,
    catalog,
    subscriptions,
    alerts,
    notebook,
    notebooks,
    dashboards,
    digests,
    reviews,
    automl,
    ml,
    chemistry,
    sphingolipid,
    phenotypes,
    graph,
    structures,
    pockets,
    docking,
    poses,
    connectivity,
    single_cell,
    spectral,
    crispr,
    multi_omics,
    multi_omics_viz,
    variants,
    companies,
    extraction,
    sync,
    admet,
    qsar,
    biomarker,
    viz3d,
    compound_target,
    explorer,
    pathway_maps,
    multi_omics_viz,
)
from amprenta_rag.config import get_config

# Get configuration
cfg = get_config()

# Create FastAPI app
app = FastAPI(
    title="Amprenta Multi-Omics Platform API",
    description="REST API for the multi-omics research platform",
    version="1.0.0",
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=cfg.server.cors_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(programs.router, prefix="/api/v1/programs", tags=["Programs"])
app.include_router(experiments.router, prefix="/api/v1/experiments", tags=["Experiments"])
app.include_router(datasets.router, prefix="/api/v1/datasets", tags=["Datasets"])
app.include_router(features.router, prefix="/api/v1/features", tags=["Features"])
app.include_router(signatures.router, prefix="/api/v1/signatures", tags=["Signatures"])
app.include_router(compounds.router, prefix="/api/v1/compounds", tags=["Compounds"])
app.include_router(screening.router, prefix="/api/v1/screening", tags=["Screening"])
app.include_router(sar.router, prefix="/api/v1/sar", tags=["SAR"])
app.include_router(
    explainability.router, prefix="/api/v1", tags=["explainability"]
)
app.include_router(analysis.router, prefix="/api", tags=["Analysis"])
app.include_router(batch.router, prefix="/api", tags=["Analysis"])
app.include_router(sphingolipid.router, prefix="/api", tags=["Analysis"])
app.include_router(reports.router, prefix="/api/v1", tags=["Reports"])
app.include_router(quality.router, prefix="/api/v1", tags=["Quality"])
app.include_router(protocols.router, prefix="/api/v1", tags=["Protocols"])
app.include_router(hts.router, prefix="/api/v1", tags=["HTS"])
app.include_router(pathways.router, prefix="/api/v1", tags=["Pathways"])
app.include_router(moa.router, prefix="/api/v1", tags=["MOA"])
app.include_router(chemistry.router, prefix="/api/v1", tags=["Chemistry"])
app.include_router(catalog.router, prefix="/api/v1", tags=["Catalog"])
app.include_router(subscriptions.router, prefix="/api/v1", tags=["Subscriptions"])
# Note: `amprenta_rag.api.routers.alerts` contains TWO routers:
# - `alerts.router` (legacy notification alerts) served at `/api/v1/alerts/*`
# - `alerts.structural_router` (structural alert checking: PAINS/Brenk/Lilly) served at `/api/alerts/*`
app.include_router(alerts.router, prefix="/api/v1", tags=["Alerts"])
app.include_router(alerts.structural_router, prefix="/api")
app.include_router(notebook.router, prefix="/api/notebook", tags=["Notebook"])
app.include_router(notebooks.router, prefix="/api")
app.include_router(dashboards.router, prefix="/api")
app.include_router(digests.router, prefix="/api")
app.include_router(reviews.router, prefix="/api")
app.include_router(automl.router, prefix="/api")
app.include_router(ml.router, prefix="/api")

app.include_router(phenotypes.router, prefix="/api")
app.include_router(graph.router, prefix="/api")
app.include_router(structures.router, prefix="/api")
app.include_router(pockets.router, prefix="/api")
app.include_router(docking.router, prefix="/api")
app.include_router(poses.router, prefix="/api")
app.include_router(connectivity.router, prefix="/api")
app.include_router(single_cell.router, prefix="/api")
app.include_router(spectral.router, prefix="/api")
app.include_router(crispr.router, prefix="/api")
app.include_router(multi_omics.router, prefix="/api")
app.include_router(multi_omics_viz.router, prefix="/api")
app.include_router(variants.router, prefix="/api")
app.include_router(companies.router, prefix="/api")
app.include_router(extraction.router, prefix="/api")
app.include_router(sync.router, prefix="/api")
app.include_router(admet.router, prefix="/api")
app.include_router(qsar.router, prefix="/api")
app.include_router(biomarker.router, prefix="/api")
app.include_router(viz3d.router, prefix="/api")
app.include_router(compound_target.router, prefix="/api/network")
app.include_router(explorer.router, prefix="/api")
app.include_router(pathway_maps.router, prefix="/api")
app.include_router(multi_omics_viz.router, prefix="/api")


@app.get("/")
async def root() -> dict:
    """Root endpoint."""
    return {
        "name": "Amprenta Multi-Omics Platform API",
        "version": "1.0.0",
        "status": "operational",
    }


@app.get("/health")
async def health_check() -> dict:
    """Health check endpoint."""
    return {"status": "healthy"}

