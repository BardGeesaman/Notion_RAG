# Production Hardening Guide

This document describes the production hardening features implemented for the Amprenta RAG system.

## Overview

Production hardening ensures the system is robust, reliable, and production-ready. This includes:

- ✅ Enhanced error handling and recovery
- ✅ Configuration validation
- ✅ Health checks
- ✅ Performance monitoring
- ✅ Retry logic for external APIs
- ✅ Graceful degradation

## Features

### 1. Configuration Validation

**Purpose**: Validate all configuration before application startup.

**Usage**:

```bash
python scripts/validate_configuration.py
```

**What it checks**:
- Required API keys (OpenAI, Pinecone, Notion, Zotero)
- Database IDs configuration
- External service connectivity
- Pipeline configuration values

**Example Output**:
```
✅ All configuration checks passed!
System is ready to use.
```

**Programmatic Usage**:

```python
from amprenta_rag.utils.config_validation import validate_configuration

if validate_configuration():
    print("Configuration is valid!")
```

### 2. Health Checks

**Purpose**: Monitor system health and external service status.

**Usage**:

```bash
python scripts/health_check.py
```

**What it checks**:
- Notion API connectivity and latency
- Pinecone API connectivity and latency
- OpenAI API connectivity and latency
- Configuration health

**Example Output**:
```
✅ Notion API: HEALTHY
   Connected successfully
   Latency: 332.2ms

✅ Pinecone API: HEALTHY
   Connected successfully
   Latency: 1063.1ms
   total_vectors: 2809

✅ System is healthy!
```

**Programmatic Usage**:

```python
from amprenta_rag.utils.health_check import check_system_health

results = check_system_health()
print(f"Overall status: {results['overall_status']}")
```

### 3. Error Handling & Retry Logic

**Purpose**: Automatic retry with exponential backoff for transient failures.

**Usage**:

```python
from amprenta_rag.utils.error_handling import retry_with_backoff, RetryConfig

# Basic retry with defaults
@retry_with_backoff()
def call_external_api():
    # Your API call here
    pass

# Custom retry configuration
config = RetryConfig(
    max_attempts=5,
    initial_delay=1.0,
    max_delay=60.0,
    exponential_base=2.0,
)

@retry_with_backoff(config=config)
def important_operation():
    pass
```

**Circuit Breaker Pattern**:

```python
from amprenta_rag.utils.error_handling import CircuitBreaker

breaker = CircuitBreaker(
    failure_threshold=5,
    recovery_timeout=60.0,
)

def call_service():
    return breaker.call(actual_api_call)
```

**Graceful Degradation**:

```python
from amprenta_rag.utils.error_handling import graceful_degradation

@graceful_degradation(fallback_value="default_value")
def optional_feature():
    # If this fails, returns "default_value"
    pass
```

### 4. Performance Monitoring

**Purpose**: Track and monitor performance metrics.

**Usage**:

```python
from amprenta_rag.utils.performance import timed, timer, get_performance_metrics

# Decorator approach
@timed(operation_name="database_query", log_threshold=1.0)
def query_database():
    pass

# Context manager approach
with timer("complex_operation"):
    result = perform_operation()

# Get metrics summary
metrics = get_performance_metrics()
summary = metrics.get_summary()
print(summary)
```

**Performance Metrics**:
- Execution time tracking
- Automatic logging for slow operations
- Aggregated statistics (mean, min, max, count)

### 5. Error Tracking

**Purpose**: Track and report errors for monitoring.

**Usage**:

```python
from amprenta_rag.utils.error_handling import get_error_tracker

tracker = get_error_tracker()

try:
    risky_operation()
except Exception as e:
    tracker.record_error(
        operation="risky_operation",
        error=e,
        context={"user_id": "123", "dataset": "ST001"},
    )

# Get error summary
summary = tracker.get_error_summary()
print(f"Total errors: {summary['total_errors']}")
```

## Best Practices

### 1. Use Retry Logic for External APIs

Wrap external API calls with retry logic:

```python
from amprenta_rag.utils.error_handling import retry_with_backoff
from amprenta_rag.utils.error_handling import RetryConfig
import requests

config = RetryConfig(
    max_attempts=3,
    initial_delay=1.0,
    retryable_exceptions=(requests.exceptions.RequestException,),
)

@retry_with_backoff(config=config)
def call_notion_api(url, headers):
    response = requests.get(url, headers=headers, timeout=30)
    response.raise_for_status()
    return response.json()
```

### 2. Validate Configuration on Startup

Always validate configuration before processing:

```python
from amprenta_rag.utils.config_validation import validate_configuration

def main():
    if not validate_configuration():
        logger.error("Configuration validation failed!")
        sys.exit(1)
    
    # Continue with application logic
```

### 3. Monitor Performance

Add performance tracking to critical operations:

```python
from amprenta_rag.utils.performance import timed

@timed(operation_name="signature_scoring", log_threshold=2.0)
def score_signatures():
    # Only logs if takes > 2 seconds
    pass
```

### 4. Use Health Checks for Monitoring

Set up periodic health checks:

```bash
# Cron job example
*/5 * * * * cd /path/to/project && python scripts/health_check.py
```

### 5. Handle Errors Gracefully

Use graceful degradation for optional features:

```python
from amprenta_rag.utils.error_handling import graceful_degradation

@graceful_degradation(fallback_value=[])
def fetch_optional_data():
    # Returns [] if this fails
    return fetch_data()
```

## Monitoring Integration

### Health Check Endpoint

For production deployments, you can expose health checks via HTTP:

```python
from flask import Flask, jsonify
from amprenta_rag.utils.health_check import check_system_health

app = Flask(__name__)

@app.route('/health')
def health():
    results = check_system_health()
    status_code = 200 if results['overall_status'] == 'healthy' else 503
    return jsonify(results), status_code
```

### Error Tracking Integration

Integrate error tracking with monitoring systems:

```python
from amprenta_rag.utils.error_handling import get_error_tracker

def send_to_monitoring():
    tracker = get_error_tracker()
    summary = tracker.get_error_summary()
    
    # Send to your monitoring system (e.g., Datadog, New Relic)
    monitoring_client.record_errors(summary)
```

## Troubleshooting

### Configuration Validation Fails

1. Check that all required API keys are set in `.env` file
2. Verify database IDs are correct
3. Check network connectivity to external services
4. Review the validation report for specific failures

### Health Checks Fail

1. Check API keys are valid
2. Verify network connectivity
3. Check service status pages for external services
4. Review latency values (high latency may indicate network issues)

### Performance Issues

1. Use performance monitoring to identify slow operations
2. Check error tracking for repeated failures
3. Review circuit breaker states
4. Consider increasing retry delays for overloaded services

## Future Enhancements

- [ ] Integrate with monitoring services (Datadog, New Relic, etc.)
- [ ] Add metrics export (Prometheus format)
- [ ] Implement distributed tracing
- [ ] Add alerting based on health check failures
- [ ] Create dashboard for monitoring

## Related Documentation

- [Configuration Guide](CONFIGURATION.md)
- [API Reference](API_REFERENCE.md)
- [Usage Examples](USAGE_EXAMPLES.md)

