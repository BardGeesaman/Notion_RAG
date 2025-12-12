"""Workflow automation engine."""
from __future__ import annotations

from datetime import datetime
from typing import Callable, Dict, List, Any

from amprenta_rag.database.models import WorkflowRule, WorkflowExecution
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

TRIGGER_TYPES = [
    "experiment_created",
    "compound_registered",
    "discovery_imported",
    "sample_transferred",
]

ACTION_REGISTRY: Dict[str, Callable] = {}


def register_action(action_type: str, handler: Callable) -> None:
    """
    Register an action handler.
    
    Args:
        action_type: Type of action (e.g., "send_notification")
        handler: Function that takes (config, context, db) and returns result dict
    """
    ACTION_REGISTRY[action_type] = handler
    logger.info("[AUTOMATION] Registered action handler: %s", action_type)


def fire_trigger(trigger_type: str, context: Dict[str, Any], db) -> List[WorkflowExecution]:
    """
    Fire a workflow trigger and execute matching rules.
    
    Args:
        trigger_type: Type of trigger (must be in TRIGGER_TYPES)
        context: Context dictionary with trigger data
        db: Database session
        
    Returns:
        List of WorkflowExecution records created
    """
    if trigger_type not in TRIGGER_TYPES:
        logger.warning("[AUTOMATION] Unknown trigger type: %s", trigger_type)
        return []
    
    # Query active rules matching trigger type
    rules = db.query(WorkflowRule).filter(
        WorkflowRule.trigger_type == trigger_type,
        WorkflowRule.is_active == True
    ).all()
    
    executions = []
    
    for rule in rules:
        # Check trigger_config conditions against context
        if not _matches_trigger_config(rule.trigger_config, context):
            continue
        
        # Execute action
        try:
            execution = execute_action(rule, context, db)
            executions.append(execution)
        except Exception as e:
            logger.error("[AUTOMATION] Failed to execute rule %s: %s", rule.id, e)
            # Log failed execution
            execution = WorkflowExecution(
                rule_id=rule.id,
                trigger_context=context,
                status="failed",
                result={"error": str(e)},
                triggered_at=datetime.utcnow(),
                completed_at=datetime.utcnow(),
            )
            db.add(execution)
            db.commit()
            executions.append(execution)
    
    return executions


def _matches_trigger_config(trigger_config: Dict[str, Any] | None, context: Dict[str, Any]) -> bool:
    """
    Check if context matches trigger configuration.
    
    Args:
        trigger_config: Filter conditions from rule
        context: Trigger context data
        
    Returns:
        True if context matches conditions
    """
    if not trigger_config:
        return True  # No conditions = match all
    
    # Simple matching: check if context values match config filters
    for key, expected_value in trigger_config.items():
        if key not in context:
            return False
        context_value = context.get(key)
        
        # Handle list/set comparisons
        if isinstance(expected_value, list):
            if context_value not in expected_value:
                return False
        elif isinstance(expected_value, dict):
            # Support operators like {"operator": "contains", "value": "..."}
            operator = expected_value.get("operator", "equals")
            value = expected_value.get("value")
            
            if operator == "equals":
                if context_value != value:
                    return False
            elif operator == "contains":
                if value not in str(context_value):
                    return False
            elif operator == "in":
                if context_value not in value:
                    return False
        else:
            if context_value != expected_value:
                return False
    
    return True


def execute_action(rule: WorkflowRule, context: Dict[str, Any], db) -> WorkflowExecution:
    """
    Execute a workflow rule action.
    
    Args:
        rule: WorkflowRule to execute
        context: Trigger context data
        db: Database session
        
    Returns:
        WorkflowExecution record
    """
    handler = ACTION_REGISTRY.get(rule.action_type)
    
    if not handler:
        raise ValueError(f"No handler registered for action type: {rule.action_type}")
    
    triggered_at = datetime.utcnow()
    
    try:
        # Execute handler
        result = handler(rule.action_config or {}, context, db)
        
        execution = WorkflowExecution(
            rule_id=rule.id,
            trigger_context=context,
            status="success",
            result=result,
            triggered_at=triggered_at,
            completed_at=datetime.utcnow(),
        )
        
        logger.info("[AUTOMATION] Executed rule %s (%s)", rule.id, rule.name)
    except Exception as e:
        execution = WorkflowExecution(
            rule_id=rule.id,
            trigger_context=context,
            status="failed",
            result={"error": str(e)},
            triggered_at=triggered_at,
            completed_at=datetime.utcnow(),
        )
        logger.error("[AUTOMATION] Action execution failed: %s", e)
        raise
    
    db.add(execution)
    db.commit()
    db.refresh(execution)
    
    return execution
