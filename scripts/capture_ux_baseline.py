#!/usr/bin/env python
"""Capture UX performance baseline for API endpoints."""

import json
import time
import statistics
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any

import httpx
from amprenta_rag.tests.integration.benchmarks import BenchmarkTracker


def capture_baseline() -> Dict[str, Any]:
    """Capture performance baseline for key API endpoints."""
    
    # Initialize tracker
    tracker = BenchmarkTracker()
    
    # API base URL
    base_url = "http://localhost:8000"
    
    # Test endpoints with representative calls
    test_cases = [
        # Jobs endpoints
        ("GET", "/api/v1/jobs", "jobs_list", {}),
        ("GET", "/api/v1/jobs/types", "jobs_types", {}),
        
        # Structures endpoints  
        ("GET", "/api/v1/structures", "structures_list", {"limit": 50}),
        ("GET", "/api/v1/structures/search", "structures_search", {"query": "protein"}),
        
        # Poses endpoints (if data exists)
        ("GET", "/api/v1/poses", "poses_list", {"limit": 20}),
        
        # Docking endpoints
        ("GET", "/api/v1/docking/runs", "docking_runs_list", {"limit": 30}),
        
        # Additional high-traffic endpoints
        ("GET", "/api/v1/programs", "programs_list", {}),
        ("GET", "/api/v1/datasets", "datasets_list", {"limit": 25}),
        ("GET", "/api/v1/compounds", "compounds_list", {"limit": 25}),
    ]
    
    measurements = {}
    
    try:
        with httpx.Client(timeout=30.0) as client:
            print("🔍 Capturing performance baseline...")
            
            for method, endpoint, test_name, params in test_cases:
                print(f"   Testing {method} {endpoint}")
                
                # Run multiple iterations for statistical accuracy
                times = []
                
                for i in range(5):  # 5 iterations per endpoint
                    start = time.perf_counter()
                    
                    try:
                        if method == "GET":
                            response = client.get(f"{base_url}{endpoint}", params=params)
                        elif method == "POST":
                            response = client.post(f"{base_url}{endpoint}", json=params)
                        else:
                            continue
                        
                        elapsed_ms = (time.perf_counter() - start) * 1000
                        
                        # Only record successful responses
                        if response.status_code in [200, 201]:
                            times.append(elapsed_ms)
                            tracker.record(f"{test_name}_iter_{i}", endpoint, method, elapsed_ms)
                        else:
                            print(f"     Warning: {response.status_code} for {endpoint}")
                            
                    except Exception as e:
                        print(f"     Error: {e} for {endpoint}")
                        continue
                
                # Calculate statistics if we have measurements
                if times:
                    measurements[test_name] = {
                        "endpoint": endpoint,
                        "method": method,
                        "samples": len(times),
                        "p50": statistics.median(times),
                        "p95": statistics.quantiles(times, n=20)[18] if len(times) >= 5 else max(times),
                        "p99": max(times),
                        "mean": statistics.mean(times),
                        "min": min(times),
                        "max": max(times)
                    }
                    print(f"     ✓ p50: {measurements[test_name]['p50']:.1f}ms")
                else:
                    print(f"     ✗ No successful measurements")
    
    except Exception as e:
        print(f"❌ Error during baseline capture: {e}")
        return {}
    
    # Generate comprehensive baseline report
    baseline_report = {
        "generated_at": datetime.utcnow().isoformat(),
        "purpose": "UX Performance Baseline - Pre-Optimization",
        "total_endpoints_tested": len(test_cases),
        "successful_measurements": len(measurements),
        "measurements": measurements,
        "summary_statistics": {
            "overall_p50": statistics.median([m["p50"] for m in measurements.values()]) if measurements else 0,
            "overall_p95": statistics.median([m["p95"] for m in measurements.values()]) if measurements else 0,
            "slowest_endpoint": max(measurements.items(), key=lambda x: x[1]["p95"])[0] if measurements else None
        },
        "benchmark_tracker_report": tracker.report()
    }
    
    return baseline_report


def main():
    """Main entry point."""
    print("📊 UX Performance Baseline Capture")
    print("=" * 50)
    
    # Capture baseline
    baseline = capture_baseline()
    
    if not baseline:
        print("❌ Failed to capture baseline")
        return 1
    
    # Save to reports directory
    reports_dir = Path("reports")
    reports_dir.mkdir(exist_ok=True)
    
    baseline_file = reports_dir / "ux_performance_baseline.json"
    with open(baseline_file, "w") as f:
        json.dump(baseline, f, indent=2)
    
    print(f"\n✅ Baseline captured: {baseline_file}")
    print(f"   Endpoints tested: {baseline['total_endpoints_tested']}")
    print(f"   Successful measurements: {baseline['successful_measurements']}")
    
    if baseline['successful_measurements'] > 0:
        summary = baseline['summary_statistics']
        print(f"   Overall p50: {summary['overall_p50']:.1f}ms")
        print(f"   Overall p95: {summary['overall_p95']:.1f}ms")
        if summary['slowest_endpoint']:
            print(f"   Slowest endpoint: {summary['slowest_endpoint']}")
    
    return 0


if __name__ == "__main__":
    exit(main())
