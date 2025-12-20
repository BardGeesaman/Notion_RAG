"""Email template generators for notifications."""
from __future__ import annotations

from typing import List, Dict, Any, Optional


def get_experiment_summary_html(experiment, datasets=None) -> str:
    """
    Generate HTML email template for experiment summary.

    Args:
        experiment: Experiment model object
        datasets: Optional list of dataset objects related to the experiment

    Returns:
        HTML string with experiment details
    """
    experiment_name = getattr(experiment, "name", "Unknown Experiment")
    description = getattr(experiment, "description", "No description available")
    design_type = getattr(experiment, "design_type", "Unknown")

    # Try to get organism from experiment, or from datasets if available
    organism = getattr(experiment, "organism", None)
    if not organism and datasets:
        # Get unique organisms from datasets
        organisms = set()
        for dataset in datasets:
            dataset_org = getattr(dataset, "organism", None)
            if dataset_org:
                if isinstance(dataset_org, list):
                    organisms.update(dataset_org)
                else:
                    organisms.add(dataset_org)
        organism = ", ".join(sorted(organisms)) if organisms else "Unknown"
    else:
        organism = organism if organism else "Unknown"
        if isinstance(organism, list):
            organism = ", ".join(organism)

    dataset_count = len(datasets) if datasets else 0

    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <style>
            body {{ font-family: Arial, sans-serif; line-height: 1.6; color: #333; }}
            .container {{ max-width: 600px; margin: 0 auto; padding: 20px; }}
            .header {{ background-color: #4CAF50; color: white; padding: 20px; border-radius: 5px 5px 0 0; }}
            .content {{ background-color: #f9f9f9; padding: 20px; border-radius: 0 0 5px 5px; }}
            .field {{ margin-bottom: 15px; }}
            .label {{ font-weight: bold; color: #555; }}
            .value {{ margin-top: 5px; color: #333; }}
            .footer {{ margin-top: 20px; padding-top: 20px; border-top: 1px solid #ddd; font-size: 12px; color: #777; }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>Experiment Summary</h1>
            </div>
            <div class="content">
                <div class="field">
                    <div class="label">Experiment Name:</div>
                    <div class="value">{experiment_name}</div>
                </div>
                <div class="field">
                    <div class="label">Description:</div>
                    <div class="value">{description}</div>
                </div>
                <div class="field">
                    <div class="label">Design Type:</div>
                    <div class="value">{design_type}</div>
                </div>
                <div class="field">
                    <div class="label">Organism:</div>
                    <div class="value">{organism}</div>
                </div>
                <div class="field">
                    <div class="label">Datasets:</div>
                    <div class="value">{dataset_count} dataset(s) associated</div>
                </div>
            </div>
            <div class="footer">
                <p>This is an automated notification from Amprenta Multi-Omics Platform.</p>
            </div>
        </div>
    </body>
    </html>
    """
    return html


def get_digest_html(activities: List[Dict[str, Any]], period: str) -> str:
    """
    Generate HTML email template for activity digest.

    Args:
        activities: List of activity dicts with keys like 'type', 'name', 'count', 'url'
        period: Time period (e.g., "Daily", "Weekly")

    Returns:
        HTML string with activity summary
    """
    # Count activities by type
    activity_counts: Dict[str, int] = {}
    activity_items: Dict[str, List[Dict[str, Any]]] = {}

    for activity in activities:
        activity_type = activity.get("type", "unknown")
        activity_counts[activity_type] = activity_counts.get(activity_type, 0) + 1

        if activity_type not in activity_items:
            activity_items[activity_type] = []
        activity_items[activity_type].append(activity)

    # Build activity sections
    activity_sections = []
    for activity_type, items in activity_items.items():
        type_label = activity_type.replace("_", " ").title()
        section_items = []
        for item in items[:10]:  # Limit to 10 items per type
            name = item.get("name", "Unknown")
            url = item.get("url", "#")
            section_items.append(f'<li><a href="{url}">{name}</a></li>')

        if len(items) > 10:
            section_items.append(f'<li><em>... and {len(items) - 10} more</em></li>')

        activity_sections.append(f"""
            <div class="activity-section">
                <h3>{type_label} ({len(items)})</h3>
                <ul>
                    {''.join(section_items)}
                </ul>
            </div>
        """)

    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <style>
            body {{ font-family: Arial, sans-serif; line-height: 1.6; color: #333; }}
            .container {{ max-width: 600px; margin: 0 auto; padding: 20px; }}
            .header {{ background-color: #2196F3; color: white; padding: 20px; border-radius: 5px 5px 0 0; }}
            .content {{ background-color: #f9f9f9; padding: 20px; border-radius: 0 0 5px 5px; }}
            .activity-section {{ margin-bottom: 20px; }}
            .activity-section h3 {{ color: #2196F3; margin-bottom: 10px; }}
            .activity-section ul {{ list-style-type: none; padding-left: 0; }}
            .activity-section li {{ padding: 5px 0; border-bottom: 1px solid #eee; }}
            .activity-section a {{ color: #2196F3; text-decoration: none; }}
            .activity-section a:hover {{ text-decoration: underline; }}
            .summary {{ background-color: #e3f2fd; padding: 15px; border-radius: 5px; margin-bottom: 20px; }}
            .footer {{ margin-top: 20px; padding-top: 20px; border-top: 1px solid #ddd; font-size: 12px; color: #777; }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>{period} Activity Digest</h1>
            </div>
            <div class="content">
                <div class="summary">
                    <h2>Summary</h2>
                    <p>Total activities: {len(activities)}</p>
                    <ul>
                        {''.join([f'<li>{label}: {count}</li>' for label, count in activity_counts.items()])}
                    </ul>
                </div>
                {''.join(activity_sections)}
            </div>
            <div class="footer">
                <p>This is an automated {period.lower()} digest from Amprenta Multi-Omics Platform.</p>
            </div>
        </div>
    </body>
    </html>
    """
    return html


def get_share_email_html(
    entity_type: str,
    entity_name: str,
    shared_by: str,
    message: Optional[str] = None,
) -> str:
    """
    Generate HTML email template for sharing notification.

    Args:
        entity_type: Type of entity being shared (e.g., "Experiment", "Compound")
        entity_name: Name of the entity
        shared_by: Username of the person sharing
        message: Optional message from the sharer

    Returns:
        HTML string for share notification
    """
    entity_type_label = entity_type.replace("_", " ").title()

    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <style>
            body {{ font-family: Arial, sans-serif; line-height: 1.6; color: #333; }}
            .container {{ max-width: 600px; margin: 0 auto; padding: 20px; }}
            .header {{ background-color: #FF9800; color: white; padding: 20px; border-radius: 5px 5px 0 0; }}
            .content {{ background-color: #f9f9f9; padding: 20px; border-radius: 0 0 5px 5px; }}
            .share-info {{ background-color: #fff3cd; padding: 15px; border-radius: 5px; margin: 15px 0; }}
            .message-box {{ background-color: white; padding: 15px; border-left: 4px solid #FF9800; margin: 15px 0; }}
            .footer {{ margin-top: 20px; padding-top: 20px; border-top: 1px solid #ddd; font-size: 12px; color: #777; }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>Shared {entity_type_label}</h1>
            </div>
            <div class="content">
                <p><strong>{shared_by}</strong> has shared a {entity_type_label.lower()} with you:</p>
                <div class="share-info">
                    <h2>{entity_name}</h2>
                </div>
                {f'<div class="message-box"><p><em>{message}</em></p></div>' if message else ''}
                <p>You can view this {entity_type_label.lower()} in the Amprenta Multi-Omics Platform.</p>
            </div>
            <div class="footer">
                <p>This is an automated notification from Amprenta Multi-Omics Platform.</p>
            </div>
        </div>
    </body>
    </html>
    """
    return html

