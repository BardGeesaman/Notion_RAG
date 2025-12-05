#!/usr/bin/env python3
"""
Fix database locks and blocking connections.

This script helps resolve database connection issues that prevent migrations.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.database.base import get_engine
from sqlalchemy import text
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def check_locks():
    """Check for blocking locks and idle connections."""
    engine = get_engine()
    
    with engine.connect() as conn:
        # Check for idle in transaction
        result = conn.execute(text('''
            SELECT 
                pid,
                usename,
                application_name,
                state,
                query_start,
                state_change,
                LEFT(query, 100) as query_snippet
            FROM pg_stat_activity
            WHERE datname = current_database()
            AND state = 'idle in transaction'
            AND pid != pg_backend_pid()
            ORDER BY state_change
        '''))
        
        idle_connections = result.fetchall()
        
        if idle_connections:
            print(f"\n⚠️  Found {len(idle_connections)} idle in transaction connections:\n")
            for conn_info in idle_connections:
                pid, user, app, state, q_start, s_change, query = conn_info
                print(f"  PID {pid}:")
                print(f"    User: {user}")
                print(f"    Application: {app or 'N/A'}")
                print(f"    State: {state}")
                print(f"    Idle since: {s_change}")
                print(f"    Query: {query}...")
                print()
            return [row[0] for row in idle_connections]  # Return PIDs
        else:
            print("✅ No idle in transaction connections found")
            return []


def check_blocking():
    """Check for blocking locks."""
    engine = get_engine()
    
    with engine.connect() as conn:
        result = conn.execute(text('''
            SELECT 
                blocked_locks.pid AS blocked_pid,
                blocking_locks.pid AS blocking_pid,
                blocking_activity.query AS blocking_query
            FROM pg_catalog.pg_locks blocked_locks
            JOIN pg_catalog.pg_stat_activity blocked_activity ON blocked_activity.pid = blocked_locks.pid
            JOIN pg_catalog.pg_locks blocking_locks 
                ON blocking_locks.locktype = blocked_locks.locktype
                AND blocking_locks.database IS NOT DISTINCT FROM blocked_locks.database
                AND blocking_locks.relation IS NOT DISTINCT FROM blocked_locks.relation
                AND blocking_locks.pid != blocked_locks.pid
            JOIN pg_catalog.pg_stat_activity blocking_activity ON blocking_activity.pid = blocking_locks.pid
            WHERE NOT blocked_locks.granted
        '''))
        
        blocks = result.fetchall()
        
        if blocks:
            print(f"\n⚠️  Found {len(blocks)} blocking locks:\n")
            for block in blocks:
                blocked_pid, blocking_pid, query = block
                print(f"  Blocked PID {blocked_pid} is blocked by PID {blocking_pid}")
                print(f"    Blocking query: {query[:100] if query else 'N/A'}...\n")
            return True
        else:
            print("✅ No blocking locks found")
            return False


def terminate_connection(pid: int):
    """Terminate a database connection by PID."""
    engine = get_engine()
    
    with engine.connect() as conn:
        try:
            conn.execute(text(f"SELECT pg_terminate_backend({pid})"))
            conn.commit()
            print(f"✅ Terminated connection PID {pid}")
            return True
        except Exception as e:
            print(f"❌ Failed to terminate PID {pid}: {e}")
            return False


def main():
    """Main function."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Fix database locks")
    parser.add_argument(
        "--check",
        action="store_true",
        help="Check for locks and blocking connections",
    )
    parser.add_argument(
        "--terminate-idle",
        action="store_true",
        help="Terminate idle in transaction connections",
    )
    parser.add_argument(
        "--terminate-pid",
        type=int,
        help="Terminate a specific connection by PID",
    )
    
    args = parser.parse_args()
    
    if args.check or not any([args.terminate_idle, args.terminate_pid]):
        print("Checking database connections...")
        idle_pids = check_locks()
        has_blocks = check_blocking()
        
        if idle_pids and not args.terminate_idle:
            print(f"\n💡 Tip: Run with --terminate-idle to close these connections")
    
    if args.terminate_idle:
        print("\nTerminating idle in transaction connections...")
        idle_pids = check_locks()
        for pid in idle_pids:
            terminate_connection(pid)
    
    if args.terminate_pid:
        terminate_connection(args.terminate_pid)


if __name__ == "__main__":
    main()

