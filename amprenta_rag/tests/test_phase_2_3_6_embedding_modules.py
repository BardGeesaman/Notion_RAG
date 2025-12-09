"""Tests for Phase 2.3-2.6: Embedding Module Notion Removal."""

import inspect


class TestEmbeddingModuleImports:
    """Test all embedding modules import successfully."""
    
    def test_metabolomics_embedding_imports(self):
        """Test metabolomics embedding module imports."""
        from amprenta_rag.ingestion.metabolomics.embedding import embed_metabolomics_dataset
        assert callable(embed_metabolomics_dataset)
    
    def test_proteomics_embedding_imports(self):
        """Test proteomics embedding module imports."""
        from amprenta_rag.ingestion.proteomics.embedding import embed_proteomics_dataset
        assert callable(embed_proteomics_dataset)
    
    def test_lipidomics_embedding_imports(self):
        """Test lipidomics embedding module imports."""
        from amprenta_rag.ingestion.lipidomics.embedding import embed_lipidomics_dataset
        assert callable(embed_lipidomics_dataset)
    
    def test_transcriptomics_embedding_imports(self):
        """Test transcriptomics embedding module imports."""
        from amprenta_rag.ingestion.transcriptomics.embedding import embed_transcriptomics_dataset
        assert callable(embed_transcriptomics_dataset)


class TestNoNotionReferences:
    """Test embedding modules don't contain Notion references."""
    
    def test_no_notion_in_metabolomics_embedding(self):
        """Test no Notion references in metabolomics embedding."""
        from amprenta_rag.ingestion.metabolomics import embedding
        source = inspect.getsource(embedding)
        assert 'from amprenta_rag.ingestion.dataset_notion_utils' not in source
        assert 'update_dataset_embedding_metadata' not in source
        assert 'notion_client' not in source.lower()
    
    def test_no_notion_in_proteomics_embedding(self):
        """Test no Notion references in proteomics embedding."""
        from amprenta_rag.ingestion.proteomics import embedding
        source = inspect.getsource(embedding)
        assert 'from amprenta_rag.ingestion.dataset_notion_utils' not in source
        assert 'update_dataset_embedding_metadata' not in source
        assert 'notion_client' not in source.lower()
    
    def test_no_notion_in_lipidomics_embedding(self):
        """Test no Notion references in lipidomics embedding."""
        from amprenta_rag.ingestion.lipidomics import embedding
        source = inspect.getsource(embedding)
        assert 'from amprenta_rag.ingestion.dataset_notion_utils' not in source
        assert 'update_dataset_embedding_metadata' not in source
        assert 'notion_client' not in source.lower()
    
    def test_no_notion_in_transcriptomics_embedding(self):
        """Test no Notion references in transcriptomics embedding."""
        from amprenta_rag.ingestion.transcriptomics import embedding
        source = inspect.getsource(embedding)
        assert 'from amprenta_rag.ingestion.dataset_notion_utils' not in source
        assert 'update_dataset_embedding_metadata' not in source
        assert 'notion_client' not in source.lower()


class TestFunctionSignatures:
    """Test embedding functions have correct signatures."""
    
    def test_metabolomics_function_signature(self):
        """Test metabolomics embedding function signature."""
        from amprenta_rag.ingestion.metabolomics.embedding import embed_metabolomics_dataset
        sig = inspect.signature(embed_metabolomics_dataset)
        params = list(sig.parameters.keys())
        
        assert 'page_id' in params
        assert 'dataset_name' in params
        assert 'metabolites' in params
    
    def test_proteomics_function_signature(self):
        """Test proteomics embedding function signature."""
        from amprenta_rag.ingestion.proteomics.embedding import embed_proteomics_dataset
        sig = inspect.signature(embed_proteomics_dataset)
        params = list(sig.parameters.keys())
        
        assert 'page_id' in params
        assert 'dataset_name' in params
        assert 'proteins' in params
    
    def test_lipidomics_function_signature(self):
        """Test lipidomics embedding function signature."""
        from amprenta_rag.ingestion.lipidomics.embedding import embed_lipidomics_dataset
        sig = inspect.signature(embed_lipidomics_dataset)
        params = list(sig.parameters.keys())
        
        assert 'page_id' in params
        assert 'dataset_name' in params
        assert 'species' in params
    
    def test_transcriptomics_function_signature(self):
        """Test transcriptomics embedding function signature."""
        from amprenta_rag.ingestion.transcriptomics.embedding import embed_transcriptomics_dataset
        sig = inspect.signature(embed_transcriptomics_dataset)
        params = list(sig.parameters.keys())
        
        assert 'page_id' in params
        assert 'dataset_name' in params
        assert 'genes' in params
        assert 'df' in params
        assert 'gene_column' in params


class TestCoreFunctionality:
    """Test core embedding functionality is intact."""
    
    def test_metabolomics_has_pinecone_logic(self):
        """Test metabolomics embedding has Pinecone upsert logic."""
        from amprenta_rag.ingestion.metabolomics import embedding
        source = inspect.getsource(embedding)
        
        assert 'get_pinecone_index' in source
        assert 'upsert' in source
        assert 'chunk_text' in source
        assert 'embed_texts' in source
    
    def test_proteomics_has_pinecone_logic(self):
        """Test proteomics embedding has Pinecone upsert logic."""
        from amprenta_rag.ingestion.proteomics import embedding
        source = inspect.getsource(embedding)
        
        assert 'get_pinecone_index' in source
        assert 'upsert' in source
        assert 'chunk_text' in source
        assert 'embed_texts' in source
    
    def test_lipidomics_has_pinecone_logic(self):
        """Test lipidomics embedding has Pinecone upsert logic."""
        from amprenta_rag.ingestion.lipidomics import embedding
        source = inspect.getsource(embedding)
        
        assert 'get_pinecone_index' in source
        assert 'upsert' in source
        assert 'chunk_text' in source
        assert 'embed_texts' in source
    
    def test_transcriptomics_has_pinecone_logic(self):
        """Test transcriptomics embedding has Pinecone upsert logic."""
        from amprenta_rag.ingestion.transcriptomics import embedding
        source = inspect.getsource(embedding)
        
        assert 'get_pinecone_index' in source
        assert 'upsert' in source
        assert 'chunk_text' in source
        assert 'embed_texts' in source

