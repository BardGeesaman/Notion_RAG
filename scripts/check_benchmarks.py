#!/usr/bin/env python
"""Check benchmark results against thresholds for CI gates."""

import json
import sys
import argparse
from pathlib import Path
from typing import Dict, Any


def load_benchmark_report(report_path: str) -> Dict[str, Any]:
    """Load benchmark report from JSON file."""
    try:
        with open(report_path) as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"ERROR: Benchmark report not found: {report_path}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"ERROR: Invalid JSON in benchmark report: {e}")
        sys.exit(1)


def check_benchmarks(report_path: str, max_fail_rate: float = 0.05) -> int:
    """
    Check benchmark results against thresholds.
    
    Args:
        report_path: Path to benchmark JSON report
        max_fail_rate: Maximum allowed failure rate (default 5%)
    
    Returns:
        0 if pass, 1 if failure rate exceeds threshold
    """
    report = load_benchmark_report(report_path)
    
    # Extract summary data
    summary = report.get("summary", {})
    total = summary.get("total", 0)
    failed = summary.get("failed", 0)
    pass_rate = summary.get("pass_rate", 0)
    
    if total == 0:
        print("⚠️  No benchmark results found - tests may not have run")
        return 0
    
    fail_rate = failed / total
    
    print(f"📊 Benchmark Results Summary:")
    print(f"   Total tests: {total}")
    print(f"   Failed: {failed}")
    print(f"   Pass rate: {pass_rate}%")
    print(f"   Fail rate: {fail_rate:.1%}")
    
    if fail_rate > max_fail_rate:
        print(f"\n❌ FAIL: {fail_rate:.1%} of tests exceeded thresholds (max {max_fail_rate:.1%})")
        
        # Show slowest operations for debugging
        slowest = report.get("slowest", [])
        if slowest:
            print(f"\n🐌 Slowest operations:")
            for i, result in enumerate(slowest[:5], 1):
                print(f"   {i}. {result['test_name']} - {result['method']} {result['endpoint']}: "
                      f"{result['elapsed_ms']}ms (threshold: {result['threshold_ms']}ms)")
        
        return 1
    
    print(f"\n✅ PASS: {fail_rate:.1%} of tests exceeded thresholds (within {max_fail_rate:.1%} limit)")
    return 0


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Check integration test benchmark results against performance thresholds",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python check_benchmarks.py results.json                    # Check with default 5%% threshold
  python check_benchmarks.py results.json --max-fail 0.1    # Allow up to 10%% failures
  python check_benchmarks.py --help                          # Show this help
        """
    )
    
    parser.add_argument(
        "report_path", 
        nargs="?",
        default="benchmark_results.json",
        help="Path to benchmark JSON report (default: benchmark_results.json)"
    )
    
    parser.add_argument(
        "--max-fail",
        type=float,
        default=0.05,
        help="Maximum allowed failure rate (default: 0.05 = 5 percent)"
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 1.0.0"
    )
    
    args = parser.parse_args()
    
    if args.max_fail < 0 or args.max_fail > 1:
        print("ERROR: --max-fail must be between 0.0 and 1.0")
        sys.exit(1)
    
    sys.exit(check_benchmarks(args.report_path, args.max_fail))


if __name__ == "__main__":
    main()
