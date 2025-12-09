# Repository Optimization Workflow

## Overview

When optimizing repository implementations (GEO, PRIDE, MetaboLights, etc.), we follow recommendations provided by domain-specific AI agents (e.g., Gemini) who are acting as expert engineers in those fields.

## Workflow

1. **Receive Expert Recommendations**: Domain-specific instructions are provided by AI agents (e.g., Gemini acting as Senior Bioinformatics Engineer, Senior Proteomics Engineer, etc.)

2. **Apply Strict Protocol Rules**: Follow the recommendations exactly as specified, including:
   - Specific libraries to use (e.g., GEOparse, pandas)
   - API endpoints and protocols
   - File format priorities
   - Error handling patterns

3. **Implement and Test**: 
   - Replace existing implementations with optimized versions
   - Test with real repository data
   - Verify improvements

4. **Document Changes**: Create documentation summarizing the optimization

## Examples

### GEO Optimization
- **Source**: Gemini (Senior Bioinformatics Engineer)
- **Recommendation**: Use GEOparse library instead of custom parsing
- **Result**: Cleaner API, automatic caching, more robust

### PRIDE Optimization
- **Source**: Gemini (Senior Proteomics Engineer)
- **Recommendation**: Priority-based file selection (mzTab → MaxQuant → Excel → TSV/CSV) with pandas parsing
- **Result**: Better file format support, more reliable extraction

## Benefits

- **Expert Knowledge**: Leverages specialized domain expertise
- **Best Practices**: Follows industry-standard approaches
- **Consistency**: Uses recommended libraries and patterns
- **Reliability**: More robust implementations

## Notes

- Always verify recommendations work in practice
- Test with real repository data before deployment
- Document both the recommendations and the implementation

